import numpy as np
from pandas import DataFrame
from numba import jit

import lib_initialization as ini
import lib_subsurface as sub
from progressbar import progressbar

# Surface energy and mass budget model for ice sheets, by Dirk van As.
# The model can be run for a single location with a time series of p, T, RH, WS, SR, and LRin,
# or for a transect for which the input variables are height-depent.
# The current version lacks:
# - reduction of sub-surface density by sub-surface melt
#
# Update November 2015, by Robert S. Fausto (RSF)
# - HIRHAM5 subsurface scheme has replaced the former one. The subsurface scheme now includes retention by capilary forces and dry snow densification.
# - The subsurface subroutine is called "subsurface_hirham". See subroutine for description of parameters and variab<=s.
# - Layers mass (or water equivalent thickness) is constant in time. Each layer has a part which is snow (snowc), ice (snic) and water (slwc).
# - sum water equivalent thickness of layer n is thickness_weq(n) = snowc(n)+snic(n)+slwc(n). This number is constant in time.
# - Partitioning into snow, water and ice varies from time to time, and so does the density of the snow fraction (the variab<= rhofirn).
# - This means that the actual thickness (as opposed to water eqiv), is:
#   thickness_act(n) = snowc(n)*(rho_w/rhofirn(n)) + snic*(rho_w/c.rho_ice) + slwc
#
# Update Spring 2016 by Peter Langen, Robert Fausto
# - New percolation scheme translated from Peter Langen's work: possibility
# to choose between a standard bucket scheme (donoDarcy =1) and a bucket
# scheme that passes to the next layer only the amount that would be
# allowed by a Darcy flow (do_no_darcy = 0).
# - New superimposed ice formulation as translated from Peter Langen's
# FORTRAN version
#
# other updates 2016-2017 by Baptiste Vandecrux
# - constants and parameter defined in external file (see Input folder)
# - routines to set initial density/temperature profiles from data (see
# IniTemperatureDensity def)
# - Tdeep changed to the mean temperature of the elevation bin over the
# study period (see PrepareTransect def)
# - parameter C_k for the fraction of potential refreezing occuring (see
# refreeze def)
# - Choices between several parametrization for fresh snow density
# (including a WS depant).
# - Choice between different definition of precipitation
# - Possibility to run it as conduction model for use as in Humphrey et al. 2012
# - Lefebvre et al. 2003 for the roughness scale of snow and ice
# - Calonne 2012 for permeability
# - Runoff according to a Darcy flow through saturated snow
# - variable layer thickness. Dynamic merging/splitting of layers.
# - initial layer thickness depant on the accumulation

# Main script for running the surface - subsurface model
# Here you can define which year you want to compute, define your
# parameters, plot output variables, evaluate model performance\.
#
# Author: Baptiste Vandecrux (bav@geus.dk)
# ========================================================================


# All parameters are defined in a csv files in Input folder, however it is
# possible to overwrite the default value by defining them again in the
# "param{kk}" struct hereunder.


def HHsubsurf(df_aws: DataFrame, c: ini.Struct):

    (
        time,
        T,
        z_T,
        WS,
        z_WS,
        RH,
        z_RH,
        pres,
        SRin,
        SRout,
        LRin,
        LRout,
        snowfall,
        rainfall,
        T_rain,
        theta,
        theta_v,
        q,
        Tsurf,
        rho_snow,
        rho_atm,
        nu,
        mu,
        L,
        SHF,
        LHF,
        Re,
        theta_2m,
        q_2m,
        ws_10m,
        meltflux,
        melt_mweq,
        sublimation_mweq,
    ) = variables_preparation(df_aws, c)
    Tsurf_obs = T * np.nan

    # Converts RH from relative to water to relative to ice
    # RH_wrt_w = RH
    # RH = RHwater2ice(RH_wrt_w, T, pres)

    # Calculates specific humidity and saturation (needs RH relative to ice!)
    [RH, q] = SpecHumSat(RH, T, pres, c)

    (
        rhofirn,
        snowc,
        snic,
        slwc,
        dgrain,
        T_ice,
        grndc,
        grndd,
        compaction,
        zrfrz,
        zsupimp,
        ts,
        zrogl,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
    ) = ini.IniVar(time, c)

    rho = rhofirn.copy()
    SRnet = rhofirn[:, 0].copy()

    # [H_comp, GF, GFsurf, GFsubsurf, melt_mweq,snowbkt_out,sublimation_mweq,
    # SMB_mweq,  ~, L, LHF, meltflux, ~, rho,  runoff,SHF, SRnet,~,  T_ice,
    # grndc, grndd , pdgrain, refreezing,  theta_2m,q_2m,ws_10m, Tsurf,
    # snowthick,z_icehorizon,Re,Ri,err]= IniVar(c)

    # c = CalculateMeanAccumulation(time,snowfall, c)

    # preparing some variables
    theta = T + z_T * c.g / c.c_pd
    theta_v = theta * (1 + ((1 - c.es) / c.es) * q)

    for k in progressbar(range(len(time))):
        # Step 1/*: Update snowthickness and instrument heights
        if k == 0:
            snowthick = np.empty((len(time)), dtype="float64")
            snowthick[0] = c.snowthick_ini
            # snow to ice transition
            z_icehorizon = int(min(c.z_ice_max, np.floor(snowthick[0] / c.dz_ice)))
        else:
            snowthick[k] = snowthick[k - 1]  # will be updated at  of time loop

        # Step 2/*: shortwave radiation balance snow & ice penetration
        thickness_m = snowc[:, k - 1] * (c.rho_water / rhofirn[:, k - 1]) + snic[
            :, k - 1
        ] * (c.rho_water / c.rho_ice)
        depth_m = np.cumsum(thickness_m, 0)
        rho[:, k] = (
            (snowc[:, k - 1] + snic[:, k - 1] + slwc[:, k - 1])
            * c.rho_water
            / thickness_m
        )

        snowbkt[k] = snowbkt[k - 1]

        SRnet, T_ice, internal_melting = SRbalance(
            SRout, SRin, z_icehorizon, snowthick, T_ice, rho, k, c
        )

        # Step 5/*:  Surface temperature calculation
        k_eff = 0.021 + 2.5e-6 * rho[:, k] ** 2

        # effective conductivity by Anderson 1976, is ok at limits
        # thickness of the first layer in m weq for thermal transfer
        thick_first_lay = snic[0, k - 1] + snowc[0, k - 1]

        # starting surface temperature solving
        if (c.solve_T_surf == 0) & (~np.isnan(Tsurf_obs[k])):
            # if user asked for using measured surface temperature instead
            # of solving for it and that there is a measurement available
            iter_max_EB = 1  # we do one iteration in the solver bellow
            Tsurf[k] = Tsurf_obs[k]  # we use measured surface temp
        else:
            iter_max_EB = c.iter_max_EB  # otherwise just using standard value

        # Prepare parameters needed for SEB
        EB_prev = 1
        dTsurf = c.dTsurf_ini  # Initial surface temperature step in search for EB=0 (C)

        # surface roughness
        if snowthick[k] > c.smallno:
            # if there is snow
            if snowbkt[k] > c.smallno:
                # fresh snow
                z_0 = c.z0_fresh_snow
            else:
                # old snow from Lefebre et al (2003) JGR
                z_0 = max(
                    c.z0_old_snow,
                    c.z0_old_snow
                    + (c.z0_ice - c.z0_old_snow) * (rho[0, k] - 600) / (920 - 600),
                )

        else:
            # ice roughness length
            z_0 = c.z0_ice

        for findbalance in range(1, iter_max_EB):
            # SENSIBLE AND LATENT HEAT FLUX
            (
                L[k],
                LHF[k],
                SHF[k],
                theta_2m[k],
                q_2m[k],
                ws_10m[k],
                Re[k],
            ) = SensLatFluxes_bulk_opt(
                WS[k],
                nu[k],
                q[k],
                snowthick[k],
                Tsurf[k],
                theta[k],
                theta_v[k],
                pres[k],
                rho_atm[k],
                z_WS[k],
                z_T[k],
                z_RH[k],
                #z_0,
                np.float64(z_0),
                c,
            )

            # SURFACE ENERGY BUDGET
            meltflux[k], Tsurf[k], dTsurf, EB_prev, stop = SurfEnergyBudget(
                SRnet,
                LRin[k],
                Tsurf[k],
                k_eff,
                thick_first_lay,
                T_ice[:, k],
                T_rain[k],
                dTsurf,
                EB_prev,
                SHF[k],
                LHF[k],
                rainfall[k],
                c,
            )

            if iter_max_EB == 1:
                # if we are using surface temperature it might have been
                # modified by SurfEnergyBudget. So we assign it again.
                Tsurf[k] = Tsurf_obs[k]

            if stop:
                break

        if (iter_max_EB != 1) & (
            (findbalance == c.iter_max_EB) & (abs(meltflux[k]) >= 10 * c.EB_max)
        ):
            print("Problem closing energy budget")

        # Step 6/*:  Mass Budget
        # in mweq
        melt_mweq[k] = meltflux[k] * c.dt_obs / c.L_fus / c.rho_water
        sublimation_mweq[k] = LHF[k] * c.dt_obs / c.L_sub / c.rho_water  # in mweq
        # positive LHF -> deposition -> dH_subl positive

        # ========== Step 7/*:  Sub-surface model ====================================
        # GF[1:c.z_ice_max] = -k_eff[1:c.z_ice_max] * (T_ice[:c.z_ice_max-1,k] - T_ice[1:c.z_ice_max, k]) / c.dz_ice
        # GFsurf[k] = - k_eff[0] * (Tsurf[k]- T_ice[1,k]) / thick_first_lay
        c.rho_fresh_snow = rho_snow

        (
            snowc[:, k],
            snic[:, k],
            slwc[:, k],
            T_ice[:, k],
            zrfrz[:, k],
            rhofirn[:, k],
            zsupimp[:, k],
            dgrain[:, k],
            zrogl[k],
            Tsurf[k],
            grndc[:, k],
            grndd[:, k],
            pgrndcapc[k],
            pgrndhflx[k],
            dH_comp[k],
            snowbkt[k],
            compaction[:, k],
        ) = sub.subsurface(
            Tsurf[k],
            grndc[:, k - 1].copy(),
            grndd[:, k - 1].copy(),
            slwc[:, k - 1].copy(),
            snic[:, k - 1].copy(),
            snowc[:, k - 1].copy(),
            rhofirn[:, k - 1].copy(),
            T_ice[:, k].copy(),
            dgrain[:, k - 1].copy(),
            snowfall[k] + sublimation_mweq[k],  # net accumulation
            rainfall[k],  # rain
            melt_mweq[k],  # melt
            c.Tdeep,
            snowbkt[k - 1].copy(),
            c,
        )

        # bulk density
        rho[:, k] = (snowc[:, k] + snic[:, k]) / (
            snowc[:, k] / rhofirn[:, k] + snic[:, k] / c.rho_ice
        )
        # refreezing[:,k] = zrfrz[:, k] + zsupimp[:, k]
        z_icehorizon = int(np.floor(snowthick[k] / c.dz_ice))

        # if k> 1:
        #     SMB_mweq[k] =  snowfall[k] - runoff[k]             + rainfall[k] + sublimation_mweq[k]

        #     # Update BV2017: With the layer-conservative model, the surface height
        #     # can be calculated outside of the sub-surface scheme assuming that the
        #     # bottom of the collumn remains at constant depth
        #     # cumulative dry compaction
        #     H_comp[k] = H_comp(k-1) + dH_comp  #in real m

        if snowthick[k] < 0:
            snowthick[k] = 0
    # rainHF = c.rho_water * c.c_w[0] * rainfall / c.dt_obs * (T_rain-Tsurf)
    
    # print("Values, when k : " + str(k))
    # print(     WS[k],
    #             nu[k],
    #             q[k],
    #             snowthick[k],
    #             Tsurf[k],
    #             theta[k],
    #             theta_v[k],
    #             pres[k],
    #             rho_atm[k],
    #             z_WS[k],
    #             z_T[k],
    #             z_RH[k],
    #             #z_0,
    #             np.float64(z_0),
    #             c,
    # )
    return (
        L,
        LHF,
        SHF,
        theta_2m,
        q_2m,
        ws_10m,
        Re,
        melt_mweq,
        sublimation_mweq,
        snowc,
        snic,
        slwc,
        T_ice,
        zrfrz,
        rhofirn,
        zsupimp,
        dgrain,
        zrogl,
        Tsurf,
        grndc,
        grndd,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
        compaction,
    )


def IniRhoSnow(T, WS, c: ini.Struct):
    # IniRhoSnow: Initialization of fresh snow density using different
    # parametreization. If modified, check that the chosen parametrization is
    # sent properly to the subsurface scheme.
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================

    mean_T = np.mean(T)

    if c.rho_snow_scheme == 0:
        rho_snow = c.fresh_snow_dens  # constant in time and place

    if c.rho_snow_scheme == 1:
        rho_snow = 625.0 + 18.7 * (mean_T - c.T_0) + 0.293 * (mean_T - c.T_0) ** 2

    # surface snow density by Reeh et al. 2005
    if c.rho_snow_scheme == 2:
        rho_snow = 481.0 + 4.834 * (mean_T - c.T_0)
    # surface snow density by Kuipers-Munneke et al. 2015
    if c.rho_snow_scheme == 3:
        rho_snow = 369.58 - 0.02298 * c.elev
    # PLA RFO fresh snow density (May 2015): BV2017 using
    # elevation-based parametrization
    if c.rho_snow_scheme == 4:
        rho_snow = 350  # default value everywhere
        ind = T >= 258.16
        rho_snow[ind] = 50 + 1.7 * (T[ind] - 258.16) ** 1.5
        ind = (T >= 258.16) & (WS >= 5)
        rho_snow[ind] = (
            50
            + 1.7 * (T[ind] - 258.16) ** 1.5
            + 25
            + 250 * (1 - np.exp(-0.2 * (WS[ind] - 5)))
        )
    # surface snow density by Liston et al., 2007
    # only works above -15 degC
    return rho_snow


def variables_preparation(df_aws: DataFrame, c: ini.Struct):
    df_aws = df_aws.interpolate()
    time = df_aws.index.values

    T = df_aws.AirTemperature1C.values + 273.15
    z_T = df_aws.HeightTemperature1m.values
    RH = df_aws.RelativeHumidity1.values
    z_RH = df_aws.HeightHumidity1m.values
    WS = df_aws.WindSpeed1ms.values
    z_WS = df_aws.HeightWindSpeed1m.values
    pres = df_aws.AirPressurehPa.values
    SRin = df_aws.ShortwaveRadiationDownWm2.values
    SRout = df_aws.ShortwaveRadiationUpWm2.values
    LRin = df_aws.LongwaveRadiationDownWm2.values
    LRout = df_aws.LongwaveRadiationUpWm2.values
    snowfall = df_aws.Snowfallmweq.values
    rainfall = df_aws.Rainfallmweq.values

    # Tsurf_obs = df_aws.Tsurf_K.values
    # Initialization of freshsnow density for both precipiTtion at the surface
    # and use in the subsurface scheme
    rho_snow = IniRhoSnow(T, WS, c)

    # Converts RH from relative to water to relative to ice
    # RH_wrt_w = RH
    # RH = RHwater2ice(RH_wrt_w, T, pres)

    # Calculates specific humidity and saturation (needs RH relative to ice!)
    [RH, q] = SpecHumSat(RH, T, pres, c)

    # Calculated precipitation types and temp
    T_rain = np.maximum(273.15, T)

    # c = CalculateMeanAccumulation(time, snowfall, c)
    # Initial value for surface temperature
    Tsurf = np.empty((len(time)), dtype="float64") 
    Tsurf[-1] = np.mean(T[:24]) - 2.4
    # The 2 m air temperature and IR skin temperature are similar during peak
    # solar irradiance, with the mean difference in temperature equal to -0.32oC
    # when incoming solar radiation is greater than 600 W m2. There is a larger
    # difference between the two during the nighttime, with 2m air temperature
    # higher than skin temperature by an average of 2.4??C when incoming
    # radiation is less than 200 Wm2.
    #  https://doi.org/10.5194/tc-12-907-2018

    # preparing some variables
    theta = T + z_T * c.g / c.c_pd
    theta_v = theta * (1 + ((1 - c.es) / c.es) * q)
    rho_atm = 100 * pres / c.R_d / T  # atmospheric density
    # dynamic viscosity of air (Pa s) (Sutherlands' formula using C = 120 K)
    mu = 18.27e-6 * (291.15 + 120) / (T + 120) * (T / 291.15) ** 1.5
    nu = mu / rho_atm

    # output variables:
    L = np.empty((len(time)), dtype="float64")
    SHF = np.empty((len(time)), dtype="float64")
    LHF = np.empty((len(time)), dtype="float64")
    Re = np.empty((len(time)), dtype="float64")
    theta_2m = np.empty((len(time)), dtype="float64")
    q_2m = np.empty((len(time)), dtype="float64")
    ws_10m = np.empty((len(time)), dtype="float64")
    meltflux = np.empty((len(time)), dtype="float64")
    melt_mweq = np.empty((len(time)), dtype="float64")
    sublimation_mweq = np.empty((len(time)), dtype="float64")
    return (
        time,
        T,
        z_T,
        WS,
        z_WS,
        RH,
        z_RH,
        pres,
        SRin,
        SRout,
        LRin,
        LRout,
        snowfall,
        rainfall,
        T_rain,
        theta,
        theta_v,
        q,
        Tsurf,
        rho_snow,
        rho_atm,
        nu,
        mu,
        L,
        SHF,
        LHF,
        Re,
        theta_2m,
        q_2m,
        ws_10m,
        meltflux,
        melt_mweq,
        sublimation_mweq,
    )


# # SEB_SMB(rho_bulk[:, k-1], snowbkt[k-1],
# def SEB_SMB(rho_bulk, snowbkt, psnowc, psnic, pslwc, prhofirn, psnowbkt, c):
#     # Step 2/*: shortwave radiation balance snow & ice penetration
#     SRnet, T_ice = SRbalance (SRout, SRin, SRnet, z_icehorizon, snowthick, T_ice, rho_bulk, c)
#     # Step 5/*:  Surface temperature calculation
#     k_eff = 0.021 + 2.5e-6*rho_bulk**2
#     # effective conductivity by Anderson 1976, is ok at limits
#     # thickness of the first layer in m weq for thermal transfer
#     thick_first_lay = psnic[0] + psnowc[0]

#     L, LHF, SHF, theta_2m, q_2m, ws_10m, Re, meltflux, Tsurf = SolveSurfaceTemperature_bulk()

#     # Step 6/*:  Mass Budget
#     # in mweq
#     melt_mweq = meltflux * c.dt_obs/c.L_fus/c.rho_water
#     sublimation_mweq = LHF * c.dt_obs/c.L_sub/c.rho_water# in mweq
#     # positive LHF -> deposition -> dH_subl positive
#     return melt_mweq, sublimation_mweq, Tsurf


#     # Step 7/*:  Sub-surface model
#     GF[1:c.z_ice_max] = -k_eff[1:c.z_ice_max] * (T_ice[:c.z_ice_max-1, k] - T_ice[1:c.z_ice_max, k])/c.dz_ice
#     GFsurf[k] =-k_eff[0] * (Tsurf[k] - T_ice[1,k]) / thick_first_lay
# #         grndhflx = GFsurf[k]
#     pTsurf = Tsurf[k]
#     ptsoil_in = T_ice[:,k]
#     zsn = snowfall[k] + sublimation_mweq[k]
#     snmel = melt_mweq[k]
#     raind = rainfall[k]
#     c.rho_fresh_snow = rho_snow[k]

#     # bulk density
#     rho[:,k]= (snowc + snic) / (snowc/rhofirn + snic/c.rho_ice)
#     refreezing[:,k] = zrfrz + supimp
#     z_icehorizon = floor(snowthick[k]/c.dz_ice)
#     return RunName, c

# def SolveSurfaceTemperature_bulk():
#     # Prepare parameters needed for SEB
#     EB_prev = 1
#     dTsurf = c.dTsurf_ini # Initial surface temperature step in search for EB=0 (C)

#     # Update BV2017: z_0 is calculated once outside the SEB loop
#     if snowthick > c.smallno:
#         # if there is snow
#         if snowbkt > c.smallno:
#             # fresh snow
#             z_0 = c.z0_fresh_snow
#         else:
#             # old snow from Lefebre et al (2003) JGR
#             z_0 = max(c.z0_old_snow,
#                       c.z0_old_snow + (c.z0_ice -c.z0_old_snow)*(rho_bulk[0] - 600)/(920 - 600))

#     else:
#         # ice roughness length
#         z_0 = c.z0_ice

#     for findbalance in range(0, c.iter_max_EB):
#         # SENSIBLE AND LATENT HEAT FLUX
#         L, LHF, SHF, theta_2m, q_2m, ws_10m, Re = \
#             SensLatFluxes_bulk(WS, nu, q, snowthick, Tsurf, theta,
#                                theta_v, pres, rho_atm, z_WS, z_T,
#                                z_RH, z_0, c)

#         # SURFACE ENERGY BUDGET
#         meltflux, Tsurf, dTsurf, EB_prev, stop = \
#             SurfEnergyBudget(SRnet, LRin, Tsurf, k_eff, thick_first_lay,
#                              T_ice, T_rain, dTsurf, EB_prev, SHF,
#                              LHF, rainfall, c)
#         if stop:
#             break
#     # loop surface energy balance
#     if (iter_max_EB != 1) & ((findbalance == c.iter_max_EB) & (abs(meltflux) >= 10*c.EB_max)):
#         print('Problem closing energy budget')
#     return L, LHF, SHF, theta_2m, q_2m, ws_10m, Re, meltflux, Tsurf


def SurfEnergyBudget(
    SRnet,
    LRin,
    Tsurf,
    k_eff,
    thick_first_lay,
    T_ice,
    T_rain,
    dTsurf,
    EB_prev,
    SHF,
    LHF,
    rainfall,
    c,
):
    # meltflux[k], Tsurf[k], dTsurf, EB_prev, stop  \
    # = SurfEnergyBudget (SRnet, LRin[k], Tsurf[k], k_eff,thick_first_lay, \
    #                     T_ice[:,k], T_rain[k], dTsurf, EB_prev, SHF[k], LHF[k], rainfall[k],c)
    # SurfEnergyBudget: calculates the surface temperature (Tsurf) and meltflux
    # from the different elements of the energy balance. The surface
    # temperature is adjusted iteratively until equilibrium between fluxes is
    # found.
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================
    stop = 0

    # SURFACE ENERGY BUDGET
    meltflux = (
        SRnet[0]
        - SRnet[1]
        + LRin
        - c.em * c.sigma * Tsurf ** 4
        - (1 - c.em) * LRin
        + SHF
        + LHF
        - (k_eff[0]) * (Tsurf - T_ice[1]) / thick_first_lay
        + c.rho_water * c.c_w * rainfall * c.dev / c.dt_obs * (T_rain - c.T_0)
    )

    if (meltflux >= 0) & (Tsurf == c.T_0):
        # stop iteration for melting surface
        stop = 1
        return meltflux, Tsurf, dTsurf, EB_prev, stop

    if abs(meltflux) < c.EB_max:
        # stop iteration when energy components in balance
        stop = 1
        meltflux = 0
        return meltflux, Tsurf, dTsurf, EB_prev, stop

    if meltflux / EB_prev < 0:
        dTsurf = 0.5 * dTsurf
        # make surface temperature step smaller when it overshoots EB=0

    EB_prev = meltflux

    # Update BV
    if meltflux < 0:
        Tsurf = Tsurf - dTsurf
    else:
        Tsurf = min(c.T_0, Tsurf + dTsurf)
    return meltflux, Tsurf, dTsurf, EB_prev, stop


def SRbalance(SRout, SRin, z_icehorizon, snowthick, T_ice, rho, k, c: ini.Struct):
    # SRbalance: Calculates the amount of Shortwave Radiation that is
    # penetrating at each layer (SRnet).
    # uses it to warm each layer and eventually calculates the melt that is
    # produced by this warming
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================
    # extinction coefficient of ice 0.6 to 1.5 m-1
    # extinction coefficient of snow 4 to 40 m-1
    # Greufell and Maykut, 1977, Journal of Glacionp.logy, Vol. 18, No. 80, 1977

    # radiation absorption in snow
    # SRnet(snow_layer) = (SRin - SRout)#     *exp(-ext_snow*depth(snow_layer))
    #
    # radiation absorption in ice layers underneath the snowpack
    #  SRnet(ice_layer) = (SRin-SRout).*#         exp(-ext_snow*snowthick).*#         exp(-ext_ice*(depth(ice_layer) - snowthick))
    SRnet = np.empty_like(T_ice[:, 0])
    # radiation absorption in snow
    SRnet[:(z_icehorizon)] = (SRin[k] - SRout[k]) * np.exp(
        -c.ext_snow * np.arange(0, z_icehorizon) * c.dz_ice
    )

    # radiation absorption in underlying ice
    if z_icehorizon < c.z_ice_max:
        SRnet[z_icehorizon:] = (
            (SRin[k] - SRout[k])
            * np.exp(-c.ext_snow * snowthick[k])
            * np.exp(
                -c.ext_ice
                * (np.arange(z_icehorizon, len(SRnet)) * c.dz_ice - snowthick[k])
            )
        )

    # specific heat of ice
    # (perhaps a slight overestimation for near-melt ice temperatures (max 48 J/kg/K))
    # snow & ice temperature rise due to shortwave radiation absorption
    # Specific heat of ice (a slight overestimation for near-melt T (max 48 J kg-1 K-1))
    c_i = 152.456 + 7.122 * T_ice[:, k - 1]

    T_ice[: c.z_ice_max - 1, k] = (
        T_ice[: c.z_ice_max - 1, k - 1]
        + c.dt_obs
        / c.dev
        / rho[: c.z_ice_max - 1, k]
        / c_i[: c.z_ice_max - 1]
        * (SRnet[: c.z_ice_max - 1] - SRnet[1 : c.z_ice_max])
        / c.dz_ice
    )

    T_ice[c.z_ice_max - 1, k] = T_ice[c.z_ice_max - 1, 0]

    # finding where/how much melt occurs
    subsurfmelt = T_ice[:, k] > c.T_0
    nosubsurfmelt = T_ice[:, k] <= c.T_0
    dT_ice = T_ice[:, k] - c.T_0
    dT_ice[nosubsurfmelt] = 0

    meltflux_internal_temp = rho[:, k] * c_i * dT_ice / c.dt_obs * c.dev * c.dz_ice
    meltflux_internal = np.sum(meltflux_internal_temp[: c.z_ice_max])

    dH_melt_internal = np.zeros_like(SRnet)
    dH_melt_internal = -c_i * dT_ice / c.L_fus * c.dz_ice
    dH_melt_internal[0] = 0

    if np.sum(subsurfmelt) > 0:
        T_ice[subsurfmelt, k] = c.T_0  # removing non-freezing temperatures

    # Reduce sub-surface density due to melt? Will most likely cause model instability
    return SRnet, T_ice, meltflux_internal


# def SRbalance (SRout, SRin, psnowc, psnic, pslwc, prhofirn, psnowbkt, T_ice, rho,  k, c):
#     # SRbalance: Calculates the amount of Shortwave Radiation that is
#     # penetrating at each layer (SRnet).
#     # uses it to warm each layer and eventually calculates the melt that is
#     # produced by this warming
#     #
#     # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
#     # translated to python by Baptiste Vandecrux (bav@geus.dk)
#     # ==========================================================================
#     # extinction coefficient of ice 0.6 to 1.5 m-1
#     # extinction coefficient of snow 4 to 40 m-1
#     # Greufell and Maykut, 1977, Journal of Glaciology, Vol. 18, No. 80, 1977

#     # radiation absorption in snow
#     # SRnet(snow_layer) = (SRin - SRout)*np.exp(-ext_snow*depth(snow_layer))
#     #
#     # radiation absorption in ice layers underneath the snowpack
#     # SRnet(ice_layer) = (SRin-SRout)* np.exp(-ext_snow*snowthick)*np.exp(-ext_ice*(depth(ice_layer) - snowthick))

#     thickness_m = psnowc * (c.rho_water / prhofirn)+ psnic * (c.rho_water / c.rho_ice)
#     depth_m = np.cumsum(thickness_m, 0)
#     rho_bulk = (psnowc+psnic+pslwc)*c.rho_water/thickness_m

#     # Step 1/*: Update snowthickness
#     z_icehorizon = np.argwhere(rho_bulk>830)
#     if len(z_icehorizon)==0:
#         z_icehorizon = len(rho_bulk)
#     if z_icehorizon == 0:
#         snowthick = psnowbkt*c.rho_water/315
#     else:
#         snowthick = psnowbkt*c.rho_water/315 + depth_m[z_icehorizon]

#     #radiation absorption in snow
#     SRnet = np.empty_like(T_ice, dtype='float32')
#     SRnet[:z_icehorizon] = (SRin[k] - SRout[k]) * np.exp(-c.ext_snow * thickness_m[:z_icehorizon])

#     #radiation absorption in underlying ice
#     if z_icehorizon == 0:
#         SRnet[(z_icehorizon+1):c.z_ice_max+1] = \
#             (SRin[k]-SRout[k])* np.exp(-c.ext_snow*snowthick[k]) \
#                 * np.exp(-c.ext_ice * (np.arrange(z_icehorizon+2, c.z_ice_max+1) * c.dz_ice - snowthick[k]))
#     elif z_icehorizon < len(SRnet):
#         SRnet[(z_icehorizon+1):c.z_ice_max+1] = \
#             (SRin[k]-SRout[k])* np.exp(-c.ext_snow*snowthick[k]) \
#                 * np.exp(-c.ext_ice * (np.arrange(z_icehorizon+2, c.z_ice_max+1) * c.dz_ice - snowthick[k]))

#     # specific heat of ice
#     # (perhaps a slight overestimation for near-melt ice temperatures (max 48 J/kg/K))
#     if k==0:
#         c_i = 152.456 + 7.122 * T_ice[:,0]
#         # snow & ice temperature rise due to shortwave radiation absorption
#     else:
#         # Specific heat of ice (a slight overestimation for near-melt T (max 48 J kg-1 K-1))
#         c_i = 152.456 + 7.122 * T_ice[:, k-1]

#         T_ice[:c.z_ice_max, k] = T_ice[:c.z_ice_max, k-1] \
#  + c.dt_obs / c.dev / rho[:c.z_ice_max, k] / c_i[:c.z_ice_max] \
#                 * (SRnet[:c.z_ice_max] - SRnet[1:(c.z_ice_max+1)]) / c.dz_ice

#         # a=(SRnet[:c.z_ice_max]-SRnet(2:(c.z_ice_max+1))) / c.dz_ice
#         T_ice[c.z_ice_max+1, k] = T_ice[c.z_ice_max+1,1]

#     # finding where/how much melt occurs
#     subsurfmelt = (T_ice[:, k] > c.T_0)
#     nosubsurfmelt = (T_ice[:, k] <= c.T_0)
#     dT_ice = T_ice[:, k] - c.T_0
#     dT_ice[nosubsurfmelt] = 0

#     meltflux_internal_temp = rho[:, k] * c_i * dT_ice / c.dt_obs * c.dev * c.dz_ice
#     meltflux_internal = np.sum(meltflux_internal_temp[:c.z_ice_max])

#     dH_melt_internal = -c_i * dT_ice / c.L_fus * c.dz_ice
#     dH_melt_internal[0] = 0

#     if np.sum(subsurfmelt) > 0:
#         T_ice[subsurfmelt, k] = c.T_0 # removing non-freezing temperatures

#     # Reduce sub-surface density due to melt? Will most likely cause model instability
#     return SRout, SRnet, T_ice, meltflux_internal, dH_melt_internal


def RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c):
    u_star = c.kappa * WS / (np.log(z_WS / z_0) - psi_m2 + psi_m1)

    # rough surfaces: Smeets & Van den Broeke 2008
    Re = u_star * z_0 / nu
    z_h = z_0 * np.exp(1.5 - 0.2 * np.log(Re) - 0.11 * (np.log(Re)) ** 2)

    if z_h < 1e-6:
        z_h = 1e-6

    z_q = z_h
    return z_h, z_q, u_star, Re

def SmoothSurf_opt(
        WS: np.float64, z_0: float, psi_m1: np.float64, 
        psi_m2: np.float64, nu: np.float64, z_WS: np.float64, c: ini.Struct
    ):
    u_star = Get_u_star(c.kappa, WS, z_WS, z_0, psi_m2, psi_m1)

    Re = u_star * z_0 / nu
    if Re <= 0.135:
        ind = 0
    elif (Re > 0.135) & (Re < 2.5):
        ind = 1
    elif Re >= 2.5:
        ind = 2
    else:
        ind = float("nan")
        print("ERROR")
        print(Re)
        print(u_star)
        print(z_WS)
        print(psi_m2)
        print(psi_m1)
        print(nu)
        print(ind)

    # smooth surfaces: Andreas 1987
    z_h, z_q = Get_zh_zq(z_0, c.ch1, c.ch2, c.ch3, ind, Re, c.cq1, c.cq2, c.cq3)

    if z_h < 1e-6:
        z_h = 1e-6
    if z_q < 1e-6:
        z_q = 1e-6

    return z_h, z_q, u_star, Re

# A function called from SmoothSurf
# Returns: Computed value of u_star
@jit(nopython=True) 
def Get_u_star(kappa: float, WS: np.float64, z_WS: np.float64, z_0: float, psi_m2: float, psi_m1: float):
    return kappa * WS / (np.log(z_WS / z_0) - psi_m2 + psi_m1)
     
# A function called from SmoothSurf
# Returns: Computed value of z_h and z_q
@jit(nopython=True)
def Get_zh_zq(z_0, ch1, ch2, ch3, ind, Re, cq1, cq2, cq3):
    z_h = z_0 * np.exp(
        ch1[ind] + ch2[ind] * np.log(Re) + ch3[ind] * (np.log(Re)) ** 2
    )
    z_q = z_0 * np.exp(
        cq1[ind] + cq2[ind] * np.log(Re) + cq3[ind] * (np.log(Re)) ** 2
    )
    return z_h, z_q


# (WS, nu, q, snowthick, Tsurf,theta, theta_v , pres,rho_atm,  z_WS, z_T, z_RH, z_0, c) = \
#     (WS[k], nu[k], q[k], snowthick[k], \
#      Tsurf[k], theta[k], theta_v[k], pres[k], rho_atm[k], \
#          z_WS[k], z_T[k], z_RH[k], z_0, c)


#@profile
def SensLatFluxes_bulk_opt(
    WS: np.float64, nu: np.float64, q: np.float64, snowthick: np.float64, 
    Tsurf: np.float64, theta: np.float64, theta_v: np.float64, pres: np.float64,
    rho_atm: np.float64, z_WS: np.float64, z_T: np.float64, z_RH: np.float64, z_0: np.float64, c: ini.Struct
):
    # SensLatFluxes: Calculates the Sensible Heat Flux (SHF), Latent Heat Fluxes
    # (LHF) and Monin-Obhukov length (L). Corrects for atmospheric stability.
    #
    # see Van As et al., 2005, The np.summer Surface Energy Balance of the High
    # Antarctic Plateau, Boundary-Layer Meteoronp.logy.
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================
   
    psi_m1 = np.float64(0)
    psi_m2 = np.float64(0)

    # will be updated later
    theta_2m = theta
    q_2m = q
    ws_10m = WS

    if WS > c.WS_lim:
        # Roughness length scales for snow or ice - initial guess
        if WS < c.smallno:
            z_h = 1e-10
            z_q = 1e-10
        else:
            if snowthick > 0:
                z_h, z_q, u_star, Re = SmoothSurf_opt(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
            else:
                z_h, z_q, u_star, Re = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)

        es_ice_surf = GetEsIceSurf(Tsurf, c.T_0, c.es_0)        
        
        q_surf = c.es * es_ice_surf / (pres - (1 - c.es) * es_ice_surf)
        L = 10e4

        if (theta >= Tsurf) & (WS >= c.WS_lim):  # stable stratification
            # correction from Holtslag, A. A. M. and De Bruin, H. A. R.: 1988, ‘Applied Modelling of the Night-Time
            # Surface Energy Balance over Land’, J. Appl. Meteorol. 27, 689–704.
            for i in range(0, c.iter_max_flux):
                psi_m1 = GetPsi_m1_stable(c.aa, c.bb, c.cc, c.dd, z_0, L)
                psi_m2 = GetPsi_m2_stable(c.aa, c.bb, c.cc, c.dd, L, z_WS)
                psi_h1 = GetPsi_h1_stable(c.aa, c.bb, c.cc, c.dd, L, z_h)
                psi_h2 = GetPsi_h2_stable(c.aa, c.bb, c.cc, c.dd, L, z_T)
                psi_q = GetPsi_q_stable(c.aa, c.bb, c.cc, c.dd, L, z_q)
                psi_q2 = GetPsi_q2_stable(c.aa, c.bb, c.cc, c.dd, L, z_RH)

                # Yttre if-sats ser inte ut att ha påverkan, men kanske ändå har det?
                if WS < c.smallno:
                    z_h = 1e-10
                    z_q = 1e-10
                else:
                    if snowthick > 0:
                        z_h, z_q, u_star, Re = SmoothSurf_opt(
                        WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )
                    else:
                        z_h, z_q, u_star, Re = RoughSurf(
                        WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )

                th_star, q_star = Get_th_star_q_star(
                    c.kappa, theta,
                    Tsurf, z_T,
                    z_h, psi_h2, psi_h1,
                    q, q_surf, z_RH,
                    z_q, psi_q2, psi_q
                )
                SHF, LHF = GetSHF_LHF(rho_atm, u_star, th_star, q_star, c.c_pd, c.L_sub)

                L_prev = L
                # L = u_star**2 * theta_v  / ( 3.9280 * th_star*(1 + 0.6077*q_star))

                L = Get_L(u_star, theta, c.es, q, c.g, c.kappa, th_star, q_star)

                if (L < c.smallno) | (abs((L_prev - L)) < c.L_dif):
                    # calculating 2m temperature, humidity and wind speed
                    theta_2m = Tsurf + th_star / c.kappa * (
                        np.log(2 / z_h) - psi_h2 + psi_h1
                    )
                    q_2m = q_surf + q_star / c.kappa * (
                        np.log(2 / z_q) - psi_q2 + psi_q
                    )
                    ws_10m = u_star / c.kappa * (np.log(10 / z_0) - psi_m2 + psi_m1)
                    break

        if (theta < Tsurf) & (WS >= c.WS_lim):  # unstable stratification
            # correction defs as in
            # Dyer, A. J.: 1974, ‘A Review of Flux-Profile Relationships’, Boundary-Layer Meteorol. 7, 363– 372.
            # Paulson, C. A.: 1970, ‘The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer’, J. Appl. Meteorol. 9, 857–861.

            for i in range(0, c.iter_max_flux):
                x1, x2, y1, y2, yq, yq2 = Compute_x_y_const(c.gamma, z_0, z_WS, z_h, z_T, z_q, z_RH, L)
                psi_m1, psi_m2, psi_h1, psi_h2, psi_q, psi_q2 = GetPsi_unstable(x1, x2, y1, y2, yq, yq2)

                if WS < c.smallno:
                    z_h = 1e-10
                    z_q = 1e-10
                else:
                    if snowthick > 0:
                        [z_h, z_q, u_star, Re] = SmoothSurf_opt(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )
                    else:
                        [z_h, z_q, u_star, Re] = RoughSurf(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )

                th_star, q_star = Get_th_star_q_star(
                    c.kappa, theta,
                    Tsurf, z_T,
                    z_h, psi_h2, psi_h1,
                    q, q_surf, z_RH,
                    z_q, psi_q2, psi_q
                )
                SHF, LHF = GetSHF_LHF(rho_atm, u_star, th_star, q_star, c.c_pd, c.L_sub)

                L_prev = L

                L = Get_L(u_star, theta, c.es, q, c.g, c.kappa, th_star, q_star)
               
                if abs((L_prev - L)) < c.L_dif:
                    # calculating 2m temperature, humidity and wind speed
                    theta_2m = Tsurf + th_star / c.kappa * (
                        np.log(2 / z_h) - psi_h2 + psi_h1
                    )
                    q_2m = q_surf + q_star / c.kappa * (
                        np.log(2 / z_q) - psi_q2 + psi_q
                    )
                    ws_10m = u_star / c.kappa * (np.log(10 / z_0) - psi_m2 + psi_m1)
                    break

    else:
        # threshold in windspeed ensuring the stability of the SHF/THF
        # caluclation. For low wind speeds those fluxes are anyway very small.
        Re = 0
        u_star = 0
        th_star = -999
        q_star = -999
        L = -999
        SHF = 0
        LHF = 0
        z_h = 1e-10
        z_q = 1e-10
        psi_m1 = 0
        psi_m2 = -999
        psi_h1 = 0
        psi_h2 = -999
        psi_q = 0
        psi_q2 = -999

        # calculating 2m temperature, humidity and wind speed
        theta_2m = theta
        q_2m = q
        ws_10m = WS

    return L, LHF, SHF, theta_2m, q_2m, ws_10m, Re

def SpecHumSat(RH, T, pres, c: ini.Struct):
    # SpecHumSat
    # - calculates saturation vapour pressure of water (es_wtr) ice (es_ice)
    # given the temperature T
    # - corrects RH for supersaturation
    # - calculates the specific humidity at saturation (qsat) for sub- and supra-
    # freezing conditions
    # - calculates the specific humidity as def of RH and qsat
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================

    # SPECIFIC HUMIDITY & SATURATION
    # saturation vapour pressure above 0 C (hPa)
    es_wtr = 10 ** (
        -7.90298 * (c.T_100 / T - 1)
        + 5.02808 * np.log10(c.T_100 / T)
        - 1.3816e-7 * (10 ** (11.344 * (1.0 - T / c.T_100)) - 1.0)
        + 8.1328e-3 * (10 ** (-3.49149 * (c.T_100 / T - 1)) - 1.0)
        + np.log10(c.es_100)
    )

    es_ice = 10 ** (
        -9.09718 * (c.T_0 / T - 1.0)
        - 3.56654 * np.log10(c.T_0 / T)
        + 0.876793 * (1.0 - T / c.T_0)
        + np.log10(c.es_0)
    )  # saturation vapour pressure below 0 C (hPa)

    #Fredrika: why can't I call this one here: es_ice = GetEsIceSurf(T, c.T_0, c.es_0) ?

    q_sat = (
        c.es * es_wtr / (pres - (1 - c.es) * es_wtr)
    )  # specific humidity at saturation (incorrect below melting point)

    freezing = (
        T < c.T_0
    )  # replacing saturation specific humidity values below melting point
    if np.sum(freezing) > 0:
        q_sat[freezing] = (
            c.es * es_ice[freezing] / (pres[freezing] - (1 - c.es) * es_ice[freezing])
        )

    # supersaturated = find(RH > 100)# replacing values of supersaturation by saturation
    # if np.sum(supersaturated) > 0
    #     RH(supersaturated) = 100
    #
    q = RH * q_sat / 100  # specific humidity in kg/kg
    return RH, q


# A function computing th_star and q_star, called from SensLatFluxes_bulk
# Returns: th_star and q_star
@jit(nopython=True)
def Get_th_star_q_star(kappa, theta, Tsurf, z_T, z_h, psi_h2, psi_h1, q, q_surf, z_RH, z_q, psi_q2, psi_q):
    th_star = (kappa * (theta - Tsurf) /
               (np.log(z_T / z_h) - psi_h2 + psi_h1))
    q_star = (kappa * (q - q_surf) /
              (np.log(z_RH / z_q) - psi_q2 + psi_q))
    return th_star, q_star

# A function computing LHF and SHF, called from SensLatFluxes_bulk
# Returns: SHF, LHF - sensible and latent heat fluxes
@jit(nopython=True)
def GetSHF_LHF(rho_atm, u_star, th_star, q_star, c_pd, L_sub):
    SHF = rho_atm * c_pd * u_star * th_star
    LHF = rho_atm * L_sub * u_star * q_star
    return SHF, LHF

# A function computing L, called from SensLatFluxes_bulk
# Parameters: Gets es, g and kappa from Struct c, u_star, theta, q, th_star, q_star
# Returns L
@jit(nopython=True)
def Get_L(u_star, theta, es, q, g, kappa, th_star, q_star):
    return (
        u_star ** 2
        * theta
        * (1 + ((1 - es) / es) * q)
        / (g * kappa * th_star * (1 + ((1 - es) / es) * q_star))
    )


# A function computing saturation vapour pressure in a faster way. 
# Parameters: Tsurf, c.T_0 as cT_0, c.es_0 as es_0
# Returns: value for es_ice_surf
@jit(nopython=True)
def GetEsIceSurf(Tsurf, cT_0, es_0):
    return (10 ** (
        -9.09718 * (cT_0 / Tsurf - 1.0)
            - 3.56654 * np.log10(cT_0 / Tsurf)
            + 0.876793 * (1.0 - Tsurf / cT_0)
            + np.log10(es_0)
        )) 

# Several functions computing values for psi in a faster way. 
# Done in separate functions to maintain correct results.
# Parameters: Values from Struct c (c.aa, c.bb, c.cc, c.dd) and z_0, L, z_WS, z_h, z_T, z_q, z_RH
# Returns: calculated values for psi_m1, psi_m2, psi_h1, psi_h2, psi_q, psi_q2
@jit(nopython=True)
def GetPsi_m1_stable(aa, bb, cc, dd, z_0, L):
    return np.float64(-(
        aa * z_0 / L
        + bb * (z_0 / L - cc / dd) * np.exp(-dd * z_0 / L)
        + bb * cc / dd
        ))

@jit(nopython=True)
def GetPsi_m2_stable(aa, bb, cc, dd, L, z_WS):
    return np.float64(-(
        aa * z_WS / L
        + bb * (z_WS / L - cc / dd) * np.exp(-dd * z_WS / L)
        + bb * cc / dd
        ))

@jit(nopython=True)
def GetPsi_h1_stable(aa, bb, cc, dd, L, z_h):
    return np.float64(-(
        aa * z_h / L
        + bb * (z_h / L - cc / dd) * np.exp(-dd * z_h / L)
        + bb * cc / dd
        ))

@jit(nopython=True)
def GetPsi_h2_stable(aa, bb, cc, dd, L , z_T) :
    return np.float64(-(
        aa * z_T / L
        + bb * (z_T / L - cc / dd) * np.exp(-dd * z_T / L)
        + bb * cc / dd
        ))

@jit(nopython=True)
def GetPsi_q_stable(aa, bb, cc, dd, L, z_q):
    return np.float64(-(
        aa * z_q / L
        + bb * (z_q / L - cc / dd) * np.exp(-dd * z_q / L)
        + bb * cc / dd
        ))

@jit(nopython=True)
def GetPsi_q2_stable(aa, bb, cc, dd, L, z_RH):
    return np.float64(-(
        aa* z_RH / L
        + bb * (z_RH / L - cc / dd) * np.exp(-dd * z_RH / L)
        + bb * cc / dd
        ))
  

# A function computing x and y parameters in a faster way.
# Parameters: gamma from Struct c (c.gamma), z_0, z_WS, z_h, z_T, z_q, z_RH and L
# Returns: the calculated parameters x1, x2, y1, y2, yq, yq2
@jit(nopython=True)
def Compute_x_y_const(gamma, z_0, z_WS, z_h, z_T, z_q, z_RH, L):
    x1 = (1 - gamma * z_0 / L) ** 0.25
    x2 = (1 - gamma * z_WS / L) ** 0.25
    y1 = (1 - gamma * z_h / L) ** 0.5
    y2 = (1 - gamma * z_T / L) ** 0.5
    yq = (1 - gamma * z_q / L) ** 0.5
    yq2 = (1 - gamma * z_RH / L) ** 0.5
    return x1, x2, y1, y2, yq, yq2

# A function updating the psi values, in a fast way.
# Parameters: x1, x2, y1, y2, yq, yq2 from SensLatFluxes_bulk
# Returns: the updated psi-values
@jit(nopython=True)
def GetPsi_unstable(x1, x2, y1, y2, yq, yq2):
    psi_m1 = np.float64(
        np.log(((1 + x1) / 2) ** 2 * (1 + x1 ** 2) / 2)
        - 2 * np.arctan(x1) + np.pi / 2
            )
    psi_m2 = np.float64(
        np.log(((1 + x2) / 2) ** 2 * (1 + x2 ** 2) / 2)
        - 2 * np.arctan(x2)
        + np.pi / 2
        )
    psi_h1 = np.float64(np.log(((1 + y1) / 2) ** 2))
    psi_h2 = np.float64(np.log(((1 + y2) / 2) ** 2))
    psi_q = np.float64(np.log(((1 + yq) / 2) ** 2))
    psi_q2 = np.float64(np.log(((1 + yq2) / 2) ** 2))
    return psi_m1, psi_m2, psi_h1, psi_h2, psi_q, psi_q2



# @jit
# def Calculate2mParam(L_prev, L, L_dif, Tsurf, th_star, kappa, z_h, z_q, z_0, psi_h1, 
#                      psi_h2, q_surf, q_star, psi_q, psi_q2, u_star, psi_m1, psi_m2):
#     # calculating 2m temperature, humidity and wind speed
#     theta_2m = Tsurf + th_star / kappa * (
#         np.log(2 / z_h) - psi_h2 + psi_h1
#     )
#     q_2m = q_surf + q_star / kappa * (
#         np.log(2 / z_q) - psi_q2 + psi_q
#     )
#     ws_10m = u_star / kappa * (np.log(10 / z_0) - psi_m2 + psi_m1)    
    
#     return theta_2m, q_2m, ws_10m



# Old functions, that have been optimized, used for comparision. To be deleted:

def SmoothSurf_old(WS, z_0, psi_m1, psi_m2, nu, z_WS, c):
    u_star = c.kappa * WS / (np.log(z_WS / z_0) - psi_m2 + psi_m1)

    Re = u_star * z_0 / nu
    if Re <= 0.135:
        ind = 0
    elif (Re > 0.135) & (Re < 2.5):
        ind = 1
    elif Re >= 2.5:
        ind = 2
    else:
        print("ERROR")
        print(Re)
        print(u_star)
        print(z_WS)
        print(psi_m2)
        print(psi_m1)
        print(nu)
        print(ind)

    # smooth surfaces: Andreas 1987
    z_h = z_0 * np.exp(
        c.ch1[ind] + c.ch2[ind] * np.log(Re) + c.ch3[ind] * (np.log(Re)) ** 2
    )
    z_q = z_0 * np.exp(
        c.cq1[ind] + c.cq2[ind] * np.log(Re) + c.cq3[ind] * (np.log(Re)) ** 2
    )

    if z_h < 1e-6:
        z_h = 1e-6
    if z_q < 1e-6:
        z_q = 1e-6

    return z_h, z_q, u_star, Re 

#@profile
def SensLatFluxes_bulk_old(
    WS: np.float64, nu: np.float64, q: np.float64, snowthick: np.float64, 
    Tsurf: np.float64, theta: np.float64, theta_v: np.float64, pres: np.float64,
    rho_atm: np.float64, z_WS: np.float64, z_T: np.float64, z_RH: np.float64, z_0: float, c: ini.Struct
):
    # SensLatFluxes: Calculates the Sensible Heat Flux (SHF), Latent Heat Fluxes
    # (LHF) and Monin-Obhukov length (L). Corrects for atmospheric stability.
    #
    # see Van As et al., 2005, The np.summer Surface Energy Balance of the High
    # Antarctic Plateau, Boundary-Layer Meteoronp.logy.
    #
    # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
    # translated to python by Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================
   
    psi_m1 = 0
    psi_m2 = 0

    # will be updated later
    theta_2m = theta
    q_2m = q
    ws_10m = WS

    if WS > c.WS_lim:
        # Roughness length scales for snow or ice - initial guess
        if WS < c.smallno:
            z_h = 1e-10
            z_q = 1e-10
        else:
            if snowthick > 0:
                z_h, z_q, u_star, Re = SmoothSurf_old(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
            else:
                z_h, z_q, u_star, Re = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)

        es_ice_surf = 10 ** (
            -9.09718 * (c.T_0 / Tsurf - 1.0)
            - 3.56654 * np.log10(c.T_0 / Tsurf)
            + 0.876793 * (1.0 - Tsurf / c.T_0)
            + np.log10(c.es_0)
        )
                
        q_surf = c.es * es_ice_surf / (pres - (1 - c.es) * es_ice_surf)
        L = 10e4

        if (theta >= Tsurf) & (WS >= c.WS_lim):  # stable stratification
            # correction from Holtslag, A. A. M. and De Bruin, H. A. R.: 1988, ‘Applied Modelling of the Night-Time
            # Surface Energy Balance over Land’, J. Appl. Meteorol. 27, 689–704.
            for i in range(0, c.iter_max_flux):
                psi_m1 = -(
                    c.aa * z_0 / L
                    + c.bb * (z_0 / L - c.cc / c.dd) * np.exp(-c.dd * z_0 / L)
                    + c.bb * c.cc / c.dd
                )
                psi_m2 = -(
                    c.aa * z_WS / L
                    + c.bb * (z_WS / L - c.cc / c.dd) * np.exp(-c.dd * z_WS / L)
                    + c.bb * c.cc / c.dd
                )
                psi_h1 = -(
                    c.aa * z_h / L
                    + c.bb * (z_h / L - c.cc / c.dd) * np.exp(-c.dd * z_h / L)
                    + c.bb * c.cc / c.dd
                )
                psi_h2 = -(
                    c.aa * z_T / L
                    + c.bb * (z_T / L - c.cc / c.dd) * np.exp(-c.dd * z_T / L)
                    + c.bb * c.cc / c.dd
                )
                psi_q = -(
                    c.aa * z_q / L
                    + c.bb * (z_q / L - c.cc / c.dd) * np.exp(-c.dd * z_q / L)
                    + c.bb * c.cc / c.dd
                )
                psi_q2 = -(
                    c.aa * z_RH / L
                    + c.bb * (z_RH / L - c.cc / c.dd) * np.exp(-c.dd * z_RH / L)
                    + c.bb * c.cc / c.dd
                )

                #Update z_h, z_q, u_star, Re with new psi values
                if WS < c.smallno:
                    z_h = 1e-10
                    z_q = 1e-10
                else:
                    if snowthick > 0:
                        z_h, z_q, u_star, Re = SmoothSurf_old(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )
                    else:
                        z_h, z_q, u_star, Re = RoughSurf(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )

                th_star = (
                    c.kappa * (theta - Tsurf) / (np.log(z_T / z_h) - psi_h2 + psi_h1)
                )
                q_star = c.kappa * (q - q_surf) / (np.log(z_RH / z_q) - psi_q2 + psi_q)
                SHF = rho_atm * c.c_pd * u_star * th_star
                LHF = rho_atm * c.L_sub * u_star * q_star

                L_prev = L
                # L = u_star**2 * theta_v  / ( 3.9280 * th_star*(1 + 0.6077*q_star))
                L = (
                    u_star ** 2
                    * theta
                    * (1 + ((1 - c.es) / c.es) * q)
                    / (c.g * c.kappa * th_star * (1 + ((1 - c.es) / c.es) * q_star))
                )

                if (L < c.smallno) | (abs((L_prev - L)) < c.L_dif):
                    # calculating 2m temperature, humidity and wind speed
                    theta_2m = Tsurf + th_star / c.kappa * (
                        np.log(2 / z_h) - psi_h2 + psi_h1
                    )
                    q_2m = q_surf + q_star / c.kappa * (
                        np.log(2 / z_q) - psi_q2 + psi_q
                    )
                    ws_10m = u_star / c.kappa * (np.log(10 / z_0) - psi_m2 + psi_m1)
                    break

        if (theta < Tsurf) & (WS >= c.WS_lim):  # unstable stratification
            # correction defs as in
            # Dyer, A. J.: 1974, ‘A Review of Flux-Profile Relationships’, Boundary-Layer Meteorol. 7, 363– 372.
            # Paulson, C. A.: 1970, ‘The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer’, J. Appl. Meteorol. 9, 857–861.

            for i in range(0, c.iter_max_flux):

                x1 = (1 - c.gamma * z_0 / L) ** 0.25
                x2 = (1 - c.gamma * z_WS / L) ** 0.25
                y1 = (1 - c.gamma * z_h / L) ** 0.5
                y2 = (1 - c.gamma * z_T / L) ** 0.5
                yq = (1 - c.gamma * z_q / L) ** 0.5
                yq2 = (1 - c.gamma * z_RH / L) ** 0.5
                
                psi_m1 = (
                    np.log(((1 + x1) / 2) ** 2 * (1 + x1 ** 2) / 2)
                    - 2 * np.arctan(x1)
                    + np.pi / 2
                )
                psi_m2 = (
                    np.log(((1 + x2) / 2) ** 2 * (1 + x2 ** 2) / 2)
                    - 2 * np.arctan(x2)
                    + np.pi / 2
                )
                psi_h1 = np.log(((1 + y1) / 2) ** 2)
                psi_h2 = np.log(((1 + y2) / 2) ** 2)
                psi_q = np.log(((1 + yq) / 2) ** 2)
                psi_q2 = np.log(((1 + yq2) / 2) ** 2)

                if WS < c.smallno:
                    z_h = 1e-10
                    z_q = 1e-10
                else:
                    if snowthick > 0:
                        [z_h, z_q, u_star, Re] = SmoothSurf_old(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )
                    else:
                        [z_h, z_q, u_star, Re] = RoughSurf(
                            WS, z_0, psi_m1, psi_m2, nu, z_WS, c
                        )

                th_star = (
                    c.kappa * (theta - Tsurf) / (np.log(z_T / z_h) - psi_h2 + psi_h1)
                )
                q_star = c.kappa * (q - q_surf) / (np.log(z_RH / z_q) - psi_q2 + psi_q)
                SHF = rho_atm * c.c_pd * u_star * th_star
                LHF = rho_atm * c.L_sub * u_star * q_star

                L_prev = L
                L = (
                    u_star ** 2
                    * theta
                    * (1 + ((1 - c.es) / c.es) * q)
                    / (c.g * c.kappa * th_star * (1 + ((1 - c.es) / c.es) * q_star))
                )

                if abs((L_prev - L)) < c.L_dif:
                    # calculating 2m temperature, humidity and wind speed
                    theta_2m = Tsurf + th_star / c.kappa * (
                        np.log(2 / z_h) - psi_h2 + psi_h1
                    )
                    q_2m = q_surf + q_star / c.kappa * (
                        np.log(2 / z_q) - psi_q2 + psi_q
                    )
                    ws_10m = u_star / c.kappa * (np.log(10 / z_0) - psi_m2 + psi_m1)
                    break

    else:
        # threshold in windspeed ensuring the stability of the SHF/THF
        # caluclation. For low wind speeds those fluxes are anyway very small.
        Re = 0
        u_star = 0
        th_star = -999
        q_star = -999
        L = -999
        SHF = 0
        LHF = 0
        z_h = 1e-10
        z_q = 1e-10
        psi_m1 = 0
        psi_m2 = -999
        psi_h1 = 0
        psi_h2 = -999
        psi_q = 0
        psi_q2 = -999

        # calculating 2m temperature, humidity and wind speed
        theta_2m = theta
        q_2m = q
        ws_10m = WS

    return L, LHF, SHF, theta_2m, q_2m, ws_10m, Re 


def SmoothSurf_opt_2(
        WS: np.float64, z_0: float, psi_m1: np.float64, 
        psi_m2: np.float64, nu: np.float64, z_WS: np.float64, c: ini.Struct
    ):
    u_star = Get_u_star(c.kappa, WS, z_WS, z_0, psi_m2, psi_m1 )

    Re = u_star * z_0 / nu
    ind = 0 if Re <= 0.135 else 1 if (Re > 0.135) and (Re < 2.5) else 2 if Re >= 2.5 else float("nan") 
    if ind == float("nan"):
        print("ERROR")
        print(Re)
        print(u_star)
        print(z_WS)
        print(psi_m2)
        print(psi_m1)
        print(nu)
        print(ind)

    # smooth surfaces: Andreas 1987
    z_h, z_q = Get_zh_zq(z_0, c.ch1, c.ch2, c.ch3, ind, Re, c.cq1, c.cq2, c.cq3)

    if z_h < 1e-6:
        z_h = 1e-6
    if z_q < 1e-6:
        z_q = 1e-6

    return z_h, z_q, u_star, Re

def SmoothSurf_opt_3(
        WS: np.float64, z_0: float, psi_m1: np.float64, 
        psi_m2: np.float64, nu: np.float64, z_WS: np.float64, c: ini.Struct
    ):
    u_star = Get_u_star(c.kappa, WS, z_WS, z_0, psi_m2, psi_m1 )

    Re = u_star * z_0 / nu
    ind = Get_ind(Re, u_star, nu, z_WS, psi_m1, psi_m2)

    # smooth surfaces: Andreas 1987
    z_h, z_q = Get_zh_zq(z_0, c.ch1, c.ch2, c.ch3, int(ind), Re, c.cq1, c.cq2, c.cq3)

    if z_h < 1e-6:
        z_h = 1e-6
    if z_q < 1e-6:
        z_q = 1e-6

    return z_h, z_q, u_star, Re

@jit(nopython=True)
def Get_ind(Re, u_star, nu, z_WS, psi_m1, psi_m2):
    if Re <= 0.135:
        ind = 0
    elif (Re > 0.135) & (Re < 2.5):
        ind = 1
    elif Re >= 2.5:
        ind = 2
    else:
        ind = np.nan
        print("ERROR")
        print(Re)
        print(u_star)
        print(z_WS)
        print(psi_m2)
        print(psi_m1)
        print(nu)
        print(ind)
    return ind