# class Dog:

#     # A simple class
#     # attribute
#     attr1 = "mammal"
#     attr2 = "dog"

#     # A sample method
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from json import load
import os

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

    def disp(self):
        for property, value in vars(self).items():
            print(property, ":", value)

def IniVar(time, c):
    df_ini = InitializationSubsurface(c)

    rhofirn = np.empty((c.num_lay, len(time)), dtype="float64")
    rho = np.empty((c.num_lay, len(time)), dtype="float64")
    snowc = np.empty((c.num_lay, len(time)), dtype="float64")
    snic = np.empty((c.num_lay, len(time)), dtype="float64")
    slwc = np.empty((c.num_lay, len(time)), dtype="float64")
    dgrain = np.empty((c.num_lay, len(time)), dtype="float64")
    tsoil = np.empty((c.num_lay, len(time)), dtype="float64")
    grndc = np.empty((c.num_lay, len(time)), dtype="float64")
    grndd = np.empty((c.num_lay, len(time)), dtype="float64")
    compaction = np.empty((c.num_lay, len(time)), dtype="float64")
    zrfrz = np.empty((c.num_lay, len(time)), dtype="float64")
    zsupimp = np.empty((c.num_lay, len(time)), dtype="float64")

    Tsurf = np.empty((len(time)), dtype="float64")
    zrogl = np.empty((len(time)), dtype="float64")
    pgrndcapc = np.empty((len(time)), dtype="float64")
    pgrndhflx = np.empty((len(time)), dtype="float64")
    dH_comp = np.empty((len(time)), dtype="float64")
    snowbkt = np.empty((len(time)), dtype="float64")
    snowthick = np.empty((len(time)), dtype="float64")

    # first time step
    rhofirn[:, -1] = df_ini.rhofirn
    snic[:, -1] = df_ini.snic
    snowc[:, -1] = df_ini.snowc
    dgrain[:, -1] = df_ini.grain_size_mm
    tsoil[:, -1] = df_ini.temp_degC
    grndc[:, -1] = tsoil[:, -1]
    snowbkt[-1] = 0
    snowthick[-1] = c.snowthick_ini
    Tsurf[-1] = 263

    return (
        rhofirn,
        rho,
        snowc,
        snic,
        slwc,
        dgrain,
        tsoil,
        grndc,
        grndd,
        compaction,
        zrfrz,
        zsupimp,
        Tsurf,
        zrogl,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
        snowthick
    )


def ImportConst(ElevGrad:float=0.1):
    # ImportConst: Reads physical, site-depant, simulation-depant and
    # user-defined parameters from a set of csv files located in the ./Input
    # folder. It stores all of them in the c structure that is then passed to
    # all defs. Structures were found to be the fastest way of
    # communicating values from one def to another.
    #
    # Author: Baptiste Vandecrux (bav@geus.dk)
    # ========================================================================


    const_phy_path = "input/constants/const_phy.csv"
    const_sim_path = "input/constants/const_sim.csv"
    const_subsurf_path = "input/constants/const_subsurf.csv"

    # Load dataframes with constants and create c containing all
    df1, df2, df3 = (pd.read_csv(const_phy_path, sep=";", header=None),
                     pd.read_csv(const_sim_path, sep=";", header=None),
                     pd.read_csv(const_subsurf_path, sep=";", header=None)
    )
    df_concat = pd.concat([df1, df2, df3])
    c = df_concat.transpose()   

    c.columns = c.iloc[0, :]
    c = c.iloc[1, :]
    c = c.apply(pd.to_numeric, errors="ignore")
    c[["ch1", "ch2", "ch3", "cq1", "cq2", "cq3"]] = c[
        ["ch1", "ch2", "ch3", "cq1", "cq2", "cq3"]
    ].apply(np.fromstring, dtype=float, sep=",")
    c = c.to_dict()
    c = Struct(**c)
    # Determine local runoff time-scale  (Zuo and Oerlemans 1996). Parameters
    # are set as in Lefebre et al (JGR, 2003) = MAR value (Fettweis pers comm)
    c.t_runoff = (c.cro_1 + c.cro_2 * np.exp(-c.cro_3 * ElevGrad)) * c.t_runoff_fact
    return c


def InitializationSubsurface(
    c, 
):  # T_obs, depth_thermistor, T_ice, time, Tsurf_ini, j, c):
    # InitializationSubsurface: Sets the initial state of the sub surface parameters:
    # - snow depth
    # - temperature profile
    # - density profile
    # Author: Baptiste Vandecrux (bav@geus.dk)
    # ==========================================================================

    # Initial density profile
    filename = c.initial_state_folder_path + c.station + "_initial_density.csv"
    if not os.path.isfile(filename):
        print('Did not find initial density profile. Using "ablation_initial_density.csv".')
        filename = c.initial_state_folder_path + "ablation_initial_density.csv"

    df_ini_dens = pd.read_csv(filename, sep=";")
    df_ini_dens.loc[df_ini_dens.density_kgm3.isnull(), "density_kgm3"] = 350

    df_ini_dens["thickness_m"] = np.insert(
        np.diff(df_ini_dens.depth_m), 0, df_ini_dens.depth_m[0]
    )
    df_ini_dens["thickness_mweq"] = (
        df_ini_dens["thickness_m"] / c.rho_water * df_ini_dens.density_kgm3
    )
    df_ini_dens["depth_weq"] = np.cumsum(df_ini_dens["thickness_mweq"])
    df_ini_dens = df_ini_dens.set_index("depth_weq")

    NumLayer = c.num_lay

    depth_mod_weq = np.insert(
        np.arange(1, NumLayer + 1) ** 4 / (NumLayer) ** 4 * c.z_max, 0, 0
    )
    # here we make sure the top layers are thick enough
    # if they are too thin we augment them to c.lim_new_lay and remove the added mass from the bottom layer so that the total depth weq is still c.z_max
    thickness_mod_weq = np.diff(depth_mod_weq)
    tmp = np.maximum(0, c.lim_new_lay - thickness_mod_weq)
    thickness_mod_weq = thickness_mod_weq + tmp - np.flip(tmp)

    depth_mod_weq = np.cumsum(np.insert(thickness_mod_weq, 0, 0))

    df_mod = pd.DataFrame(depth_mod_weq[1:], columns=["depth_weq"])
    df_mod = df_mod.set_index("depth_weq")

    df_ini_dens = pd.concat([df_ini_dens, df_mod]).sort_index()
    df_ini_dens["density_kgm3"] = (
        df_ini_dens["density_kgm3"].fillna(method="bfill").values
    )

    df_ini_dens["thickness_mweq"] = np.insert(
        np.diff(df_ini_dens.index), 0, df_ini_dens.index[0]
    )
    df_ini_dens["thickness_m"] = (
        df_ini_dens.thickness_mweq * c.rho_water / df_ini_dens.density_kgm3
    )
    df_ini_dens["depth_m_2"] = np.cumsum(df_ini_dens["thickness_m"])

    # finding within which final depth bin each layer of the merged array falls in
    df_ini_dens = df_ini_dens.loc[df_ini_dens.index <= depth_mod_weq.max(), :]
    df_ini_dens["placings"] = np.digitize(df_ini_dens.index, depth_mod_weq, right=True)

    # within each final depth bin, making the average of densities weighted by the thikcness of the layers that compose each final bin
    wm = lambda x: np.average(x, weights=df_ini_dens.loc[x.index, "thickness_m"])
    df_mod["density_kgm3"] = (
        df_ini_dens.groupby("placings")
        .agg(weighted_density=("density_kgm3", wm))
        .values
    )

    if (
        df_mod["density_kgm3"].last_valid_index()
        < df_mod["density_kgm3"].index.values[-1]
    ):
        tmp = df_mod.loc[df_mod.density_kgm3.notnull(), "density_kgm3"]
        x = np.around(np.append(tmp.index.values, [30, 70]), 4)
        y = np.around(np.append(tmp.values, [830, 830]), 4)
        fo = np.poly1d(np.polyfit(x, y, 2))
        df_mod.loc[df_mod.density_kgm3.isnull(), "density_kgm3"] = fo(
            df_mod.loc[df_mod.density_kgm3.isnull(), "density_kgm3"].index.values
        )

        df_mod["density_kgm3"] = np.minimum(
            np.maximum(300, df_mod["density_kgm3"].values), 900
        )

    # ind_last = df_mod.density_kgm3.last_valid_index()
    # df_mod.loc[ind_last:, 'density_kgm3'] = 917

    df_mod["thickness_mweq"] = np.diff(depth_mod_weq)
    df_mod["thickness_m"] = (
        df_mod["thickness_mweq"] * c.rho_water / df_mod["density_kgm3"]
    )

    df_mod["depth_m"] = np.cumsum(df_mod.thickness_m)

    df_mod["rhofirn"] = df_mod.density_kgm3
    df_mod["snowc"] = df_mod["thickness_mweq"]
    df_mod["snic"] = 0

    # Initial temperature profile
    filename = c.initial_state_folder_path + c.station + "_initial_temperature.csv"
    if not os.path.isfile(filename):
        print('Did not find initial temperature profile. Using "ablation_initial_temperature.csv".')
        filename = c.initial_state_folder_path + "ablation_initial_temperature.csv"
        
    df_ini_temp = pd.read_csv(filename, sep=";")
    df_ini_temp = df_ini_temp.loc[df_ini_temp.depth_m >= 0, :]

    if df_ini_temp.depth_m.max() < df_mod.depth_m.max():
        tmp = df_ini_temp.iloc[-1, :].copy()
        tmp.depth_m = df_mod.depth_m.max()
        df_ini_temp = pd.concat([df_ini_temp, tmp.to_frame().T], ignore_index=True)

    df_mod["temp_degC"] = np.interp(
        df_mod.depth_m, df_ini_temp.depth_m, df_ini_temp.temperature_degC
    )
    df_mod["temp_degC"] = df_mod["temp_degC"].fillna(method="bfill").values + c.T_0

    # Initial grain size
    filename = c.initial_state_folder_path + 'all_sites_initial_grain_size.csv'

    df_ini_gs = pd.read_csv(filename, sep=";")
    df_ini_gs = df_ini_gs.set_index("depth_m")

    df_mod["grain_size_mm"] = (
        df_ini_gs.groupby(
            pd.cut(df_ini_gs.index, np.insert(df_mod.depth_m.values, 0, 0))
        )
        .mean()
        .values
    )
    df_mod["grain_size_mm"] = df_mod["grain_size_mm"].interpolate().values
    df_mod["grain_size_mm"] = df_mod["grain_size_mm"].fillna(method="bfill").values

    # Initial water content
    df_mod["slwc"] = 0

    if c.verbose == 1:
        fig, ax = plt.subplots(1, 4, sharey=True)
        ax = ax.flatten()
        ax[0].step(
            df_mod.density_kgm3, -df_mod.depth_m, where="pre", label="interpolated"
        )
        ax[0].step(
            df_ini_dens.density_kgm3,
            -df_ini_dens.depth_m,
            where="pre",
            label="original",
        )
        ax[0].set_xlabel("density_kgm3")
        ax[1].step(
            df_mod.temp_degC - c.T_0, -df_mod.depth_m, where="pre", label="interpolated"
        )
        ax[1].step(df_ini_temp.temperature_degC, -df_ini_temp.depth_m, label="original")
        ax[1].set_xlabel("temp_degC")
        ax[2].step(
            df_mod.grain_size_mm, -df_mod.depth_m, where="pre", label="interpolated"
        )
        ax[2].step(
            df_ini_gs.grain_size_mm, -df_ini_gs.index, where="pre", label="original"
        )
        ax[2].set_xlabel("grain_size_mm")
        ax[2].legend()
        ax[3].step(df_mod.slwc, -df_mod.depth_m, where="pre", label="interpolated")
        ax[3].set_xlabel("slwc")

    return df_mod


# def OutputName(c)
# if c.ConductionModel == 1
#     RunName = sprintf('#s_#i_ConductionOnly',        tag, c.year)
# else
#     RunName = c.station
#     RunName = [RunName, sprintf('_#i',c.year)]
#     if c.calc_CLliq == 1
#         text_Si = 'CL'
#     else
#         text_Si = sprintf('#0.2f',c.liqmax)

#     RunName = [RunName, '_IWC_', text_Si]
#     RunName = [RunName, sprintf('_#i_layers',c.jpgrnd-1)]

#     c.OutputFolder = sprintf('#s/#s',c.OutputRoot,RunName)
#     [~,~,id] =  mkdir(c.OutputFolder)
#     count = 1
#     while ~isempty(strfind(id,'DirectoryExists'))
#        count =count+1

#        c.OutputFolder = sprintf('./Output/#s_#i',RunName,count)
#        [~,~,id] =  mkdir(c.OutputFolder)

#     if count>1
#            RunName = sprintf('#s_#i',RunName,count)

# return RunName, c

# def RenamingVariables(data_out,c)

# # time
#     year = data_out.Year
#     hour = data_out.HourOfDayUTC
#     day = data_out.DayOfYear
# # leap years
#     time = year + (day + hour/24)/365
#     leapyear = find(year/4 == floor(year/4))
#     if sum(leapyear) >0
#         time(leapyear) = year(leapyear)+(day(leapyear)+hour(leapyear)/24.)/366


# # temperature, humidity and wind speed
#     if sum(strcmp(data_out.Properties.VariableNames,'AirTemperatureC'))
#         disp('Only one level was detected on the weather station.')
#         T2 = data_out.AirTemperatureC
#         T1 = NaN(size(T2))
#         RH2 = data_out.RelativeHumidity
#         RH1 = NaN(size(T2))
#         WS2 = data_out.WindSpeedms
#         WS1 = NaN(size(T1))

#         o_T1 = NaN(size(T1))
#         o_RH1 = NaN(size(T1))
#         o_WS1 = NaN(size(T1))
#         z_T1 = NaN(size(T1))
#         z_RH1 = NaN(size(T1))
#         z_WS1 = NaN(size(T1))

# # assigning height
#         if sum(strcmp('HeightWindSpeedm',data_out.Properties.VariableNames))
#             z_WS2 = data_out.HeightWindSpeedm
#             z_T2 = data_out.HeightTemperaturem
#             z_RH2 = data_out.HeightHumiditym
#         else
#             disp('Assigning measurement height from HeightSensorBoomm_raw field.')
#             z_WS2 = data_out.HeightSensorBoomm_raw + 0.4
#             z_T2 = data_out.HeightSensorBoomm_raw - 0.12
#             z_RH2 = data_out.HeightSensorBoomm_raw - 0.12


# # assigning origin
#         if sum(strcmp('WindSpeed1ms_Origin',data_out.Properties.VariableNames))
#             o_WS2 = data_out.WindSpeed1ms_Origin
#             o_T2 = data_out.AirTemperature1C_Origin
#             o_RH2 = data_out.RelativeHumidity1_Origin
#         else
#             disp('No WindSpeed1ms_Origin field specified')
#             o_T2 = zeros(size(T1))
#             o_RH2 = zeros(size(T1))
#             o_WS2 = zeros(size(T1))


#     elseif sum(strcmp(data_out.Properties.VariableNames,'AirTemperature1C'))
#         disp('Two levels detected on the weather station')
#         T1 = data_out.AirTemperature1C
#         T2 = data_out.AirTemperature2C
#         RH1 = data_out.RelativeHumidity1
#         RH2 = data_out.RelativeHumidity2
#         WS1 = data_out.WindSpeed1ms
#         WS2 = data_out.WindSpeed2ms

#         o_T1 = data_out.AirTemperature1C_Origin
#         o_T2 = data_out.AirTemperature2C_Origin
#         o_RH1 = data_out.RelativeHumidity1_Origin
#         o_RH2 = data_out.RelativeHumidity2_Origin
#         o_WS1 = data_out.WindSpeed1ms_Origin
#         o_WS2 = data_out.WindSpeed2ms_Origin

#         z_T1 = data_out.HeightTemperature1m
#         z_T2 = data_out.HeightTemperature2m
#         z_RH1 = data_out.HeightHumidity1m
#         z_RH2 = data_out.HeightHumidity2m
#         z_WS1 = data_out.HeightWindSpeed1m
#         z_WS2 = data_out.HeightWindSpeed2m
#     else
#         error('Cannot recognize temperature field in weather data file.')

#     T1 = T1 + c.T_0# Temperature in degrees Kelvin
#     T2 = T2 + c.T_0# Temperature in degrees Kelvin

# # radiation
#     LRin = data_out.LongwaveRadiationDownWm2
#     LRout = data_out.LongwaveRadiationUpWm2

#     if sum(strcmp(data_out.Properties.VariableNames,'ShortwaveRadiationDown_CorWm2'))
#         disp('Using ShortwaveRadiation***_CorWm2 field.')
#         SRin = data_out.ShortwaveRadiationDown_CorWm2
#         SRout = data_out.ShortwaveRadiationUp_CorWm2
#     else
#         SRin = data_out.ShortwaveRadiationDownWm2
#         SRout = data_out.ShortwaveRadiationUpWm2


# # other variables
#     pres = data_out.AirPressurehPa
#     Surface_Height = data_out.SurfaceHeightm
#     Tsurf_obs = min(c.T_0, ((LRout-(1-c.em)*LRin)/c.em/c.sigma).**0.25)
#     Tsurf_obs(or(isnan(LRout),isnan(LRin))) = NaN

#     ind = strfind(data_out.Properties.VariableNames,'IceTemperature')
#     ind = find(~cellfun('isempty', ind))
#     ind2 = strfind(data_out.Properties.VariableNames,'DepthThermistor')
#     ind2 = find(~cellfun('isempty', ind2))
#     num_therm = length(ind)
#     T_ice_obs = NaN(length(T1),num_therm)
#     depth_thermistor = NaN(length(T1),num_therm)

#     if ~isempty(ind2)
#         for i = 1:length(ind)
#             T_ice_obs(:,i) = data_out.(data_out.Properties.VariableNames{ind(i)})
#             depth_thermistor(:,i) = data_out.(data_out.Properties.VariableNames{ind2(i)})


#     ind =  (LRout>316)
#     if sum(ind)>0
#         if c.verbose == 1
#         fprintf('Warning: Observed surface temperature higher than 0degC\n')

# #     before = LRout(ind)
# #     LRout(ind) = LRout(ind) - (20/15 * (T(ind)-c.T_0))
# #     figure
# #     scatter(T(ind),before,'xr')
# #     hold on
# #     scatter(T(ind),LRout(ind),'ob')
# #     leg('before correction', 'after correction')
# #     xlabel('Temperature (deg C)')
# #     ylabel('Outgoing long-wave radiation (W/m**2)')
# #     title('Observed LRout > black body at 0degC')

# return time, year, day, hour, pres,    T1, T2, z_T1, z_T2, o_T1,o_T2,     RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2,     WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,    SRin,SRout, LRin, LRout, T_ice_obs,     depth_thermistor, Surface_Height, Tsurf_obs

# def ResetTemp(depth_thermistor, LRin, LRout, T_obs, rho, T_ice, time, k, c)

# #  resets surface temperature for the conduction-only model
# if ~isnan(LRout(k)) && ~isnan(LRin(k))
#     Tsurf_reset = ((LRout(k) - (1-c.em)*LRin(k)) /(c.em*c.sigma))**(1/4)
# else
#     Tsurf_reset = NaN


# # if there is thermistor record for the first time step, then reads
# # initial subsurface conditions from AWS data
#     T_ice_reset = NaN(c.jpgrnd,1)
#     if sum(~isnan(T_obs(k,:)))>1
#         depth = depth_thermistor(k,(depth_thermistor(k,:)~=0))'
#         oldtemp = T_obs(k,(depth_thermistor(k,:)~=0))'


# # calculates the new depth scale in mweq
#     depth_weq = cumsum(c.cdel)

# # calculates the new depth scale in real m
# #     depth_act = depth_weq .*c.rho_water ./rho(:,k)
#     depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k))

# # Here we add an initial surface temperature
#     depth_act = [0 depth_act]
#     depth_weq = [0 depth_weq]
#     if depth(1) ~= 0
# #if the input density profile does not contain surface temperature,
# #then we use Tsurf_in to initiate the temperature profile
#         depth = [0 depth]
#         oldtemp = [Tsurf_reset - c.T_0 oldtemp]

# # the old scale is converted from m to mweq by interpolation
#     oldscale_weq = interp1(depth_act,depth_weq,depth)

# # we limit the new scale to the values that are given within the oldscale
#     newscale_weq = depth_weq(depth_weq<=oldscale_weq())

# # the temperature is interpolated on each stamp of the new scale
#     newtemp = interp1(oldscale_weq,oldtemp,newscale_weq)
#     newtemp = newtemp + c.T_0#going back to K

# # giving the observed temperature (on appropriate scale) as initial value
# # for the subsurface temperature profile
# # There might be a shift to apply deping on whether the first value in
# # subsurface column represent the temp at depth 0 (=Tsurf) or at depth
# # c.cdel(1). To be asked.
#     T_ice_reset(1:length(newtemp)) = newtemp

# # fills the rest of the temperature profile withan interpolation that
# # s with Tdeep at the bottom of the profile

#     if length(newtemp)<length(c.cdel)
#         d1 = length(newtemp)
#         T1 = newtemp(d1)
#         d2 = c.jpgrnd
#         T2 = c.Tdeep_AWS + c.T_0
#         ind = length(newtemp)+1:(c.jpgrnd-1)
#         T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) =             (T1-T2)/(d2-d1)**2*(d2-ind).**2 + T2
#         T_ice_reset(c.jpgrnd) = c.Tdeep_AWS+ c.T_0

# #         T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) = NaN


# # removing non-freezing temperatures (just in case)
# subsurfmelt = find(T_ice_reset(:) > c.T_0)
# if sum(subsurfmelt )>0
#     T_ice_reset(subsurfmelt) = c.T_0

# # Comes from outside the def just after it
#                     depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k))
#                     depth_act = [0 depth_act]


#                 T_ice(~isnan(T_reset),k,j) = T_reset(~isnan(T_reset))


#                     [zso_capa, zso_cond] = ice_heats (c)
#                     [grndc, grndd, ~, ~]                        = update_tempdiff_params (rho(:,k), Tdeep(j)                                            , snowc, snic, T_ice(:,k,j), zso_cond, zso_capa, c)

# return Tsurf_reset, T_ice_reset

# def RHice2water(RH_wrt_i, T, pres)
# # def converting humidity with regards to ice into humidity with
# # regards to water using GOFF-GRATCH (1945) formula.
# #
# #   RH_wrt_i is an single value or an array of relative humidity with
# #   regards to ice given in percent
# #
# #   T is the corresponding air temperature in degree Celsius
# #
# #   pres is the corresponding air pressure in hPa (not used in the current
# #   form of the def)

# T_0 = 273.15
# T_100 = 373.15
# # eps is the ratio of the molecular weights of water and dry air
# eps = 0.62198

# # Rv = 461.5
# # es_wtr = eps * exp( 1/Rv * ((L + T_0 * beta)*(1/T_0 - 1/T) - beta* np.log(T./T_0)))

# # GOFF-GRATCH 1945 equation
# es_wtr = 10.**(     -7.90298*(T_100./T - 1) + 5.02808 * np.log10(T_100./T) \# saturation vapour pressure above 0 C (hPa)
#     - 1.3816E-7 * (10.**(11.344*(1.-T./T_100))-1.)     + 8.1328E-3*(10.**(-3.49149*(T_100./T-1)) -1.) + np.log10(1013.246) )

# # WMO 2012 equation (based on Goff 1957)
# # es_wtr = 10.**(#     10.79574*(1 - T_0./T) + 5.028 * np.log10(T / T_0) \# saturation vapour pressure above 0 C (hPa)
# #     + 1.50475E-4 * (1 - 10.**(-8.2969 * (T./T_0 - 1)))#     + 0.42873E-3*(10.**(4.76955*(1 - T_0./T)) -1.) +  0.78614 + 2.0 )

# es_ice = 10.**(     -9.09718 * (T_0 ./ T - 1.) - 3.56654 * np.log10(T_0 ./ T) +     0.876793 * (1. - T ./ T_0) + np.log10(6.1071)  )# saturation vapour pressure below 0 C (hPa)

# # es_ice = 10.**(#     -9.09685 * (T_0 ./ T - 1.) - 3.56654 * np.log10(T_0 ./ T) +#     0.87682 * (1. - T ./ T_0) + 0.78614 + 2.0  )

# # q_sat_wtr = eps * es_wtr./(pres-(1-eps)*es_wtr)# specific humidity at saturation over water
# # q_sat = eps * es_ice./(pres-(1-eps)*es_ice)# specific humidity at saturation over ice
# freezing=ones(size(T))
# freezing(T>=T_0) = 0
# freezing = freezing==1

# RH_wrt_w = RH_wrt_i
# RH_wrt_w(freezing) = RH_wrt_i(freezing) .*es_ice(freezing) ./ es_wtr (freezing)# specific humidity in kg/kg
# RH_wrt_w(RH_wrt_w<0)=0
# RH_wrt_w(RH_wrt_w>100) = 100

# return RH_wrt_w

# def spechum2relhum(T, pres, q, c)

# # spechum2relhum
# #
# # Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# # translated to python by Baptiste Vandecrux (bav@geus.dk)
# #==========================================================================

# # SPECIFIC HUMIDITY & SATURATION -----------------------------------------------------------------------
#     es_wtr = 10.**(-7.90298*(c.T_100./T-1) + 5.02808 * np.log10(c.T_100./T) \# saturation vapour pressure above 0 C (hPa)
#         - 1.3816E-7 * (10.**(11.344*(1.-T./c.T_100))-1.)         + 8.1328E-3*(10.**(-3.49149*(c.T_100./T-1)) -1.) + np.log10(c.es_100))

#     es_ice = 10.**(-9.09718 * (c.T_0 ./ T - 1.) - 3.56654 * np.log10(c.T_0 ./ T) +         0.876793 * (1. - T ./ c.T_0) + np.log10(c.es_0))# saturation vapour pressure below 0 C (hPa)

#     q_sat = c.es * es_wtr./(pres-(1-c.es)*es_wtr)# specific humidity at saturation (incorrect below melting point)

#     freezing = find(T < c.T_0)# replacing saturation specific humidity values below melting point
#     RH_wrt_w = q ./ q_sat*100
#     if sum(freezing) > 0
#         q_sat(freezing) = c.es * es_ice(freezing)./(pres(freezing)-(1-c.es)*es_ice(freezing))


#     RH_wrt_i = q ./ q_sat*100
#     supersaturated = find(RH_wrt_i > 100)# replacing values of supersaturation by saturation

#     if sum(supersaturated) > 0
#         RH_wrt_i(supersaturated) = 100

# return RH_wrt_i, RH_wrt_w


# def ConvertToGivepthScale(depth_old, var_old, depth_new,opt)
# # Interpolates the depth profile of a given variable old_var
# # (temperature grain size\.) according to a given scale new_depth in real m

#     transpose_at_the_ = 0

#     if isrow(depth_old) ~= isrow(depth_new)
#         error('Old and new scale should be both rows or both columns.')
#     elseif ~isrow(depth_old)
#         transpose_at_the_ = 1
#        depth_old=depth_old'
#        var_old=var_old'
#        depth_new = depth_new'

#     var_new = NaN(size(depth_new))

#     switch opt
#         case 'linear'
# # the variable is interpolated on each stamp of the new scale
#         var_new = interp1(depth_old,var_old,depth_new,'linear','extrap')

#         case 'intensive'
# # intensive variables do not dep on the size of the system
# # f.e. if you know the density of a core section and you cut it
# # in two, the two sub-sections can be assigned the same density
# # as the original section. However we need to take into account
# # into the re-sampling the case when a new section is composed
# # of two sections in the old scale. Then the density of that
# # new section is the thickness-weighted average of the two
# # original sections.
# # example: depth_old = 1:4 density_old = 100:100:400
# #           depth_new = [0.1 0.2 0.6 1.2 3.5]
# # density_new = [100.0000  100.0000  100.0000  133.3333
# # 286.9565]

#             if depth_new()>depth_old()
#                 depth_old = [depth_old, depth_new()]
#                 var_old = [var_old, var_old()]

#             left_neighbour_in_new = depth_new
#             left_neighbour_in_new(1) = 0
#             left_neighbour_in_new(2:) = depth_new(1:-1)

#             left_neighbour_in_old =  interp1([0 depth_old],[0 depth_old],depth_new,'previous')

#             ind_type_1 = left_neighbour_in_new >= left_neighbour_in_old
#             ind_type_2 = find(left_neighbour_in_new < left_neighbour_in_old)

#             var_new(ind_type_1) = interp1([0 depth_old],                [var_old(1) var_old],                depth_new(ind_type_1),'next')

#             depth_merged = [depth_old depth_new]
#             var_merged = [var_old var_new]

#             [depth_merged, ind_sorted] = sort(depth_merged)
#             var_merged = var_merged(ind_sorted)
#             var_merged(isnan(var_merged)) = interp1([0  depth_merged(~isnan(var_merged))],                [var_merged(1) var_merged(~isnan(var_merged))],                depth_merged(isnan(var_merged)),'next')

#             thick_merged = depth_merged
#             thick_merged(2:) = depth_merged(2:) - depth_merged(1:-1)

#             for i = ind_type_2
#                 i_in_merged = discretize(depth_merged,                     [left_neighbour_in_new(i) depth_new(i)])
#                 i_in_merged(isnan(i_in_merged))= 0
#                 if i~=1
#                     i_in_merged(find(i_in_merged,1,'first')) = 0# transforming the first one into 0

#                 var_new(i) = sum(var_merged(i_in_merged==1).*thick_merged(i_in_merged==1))                     ./sum(thick_merged(i_in_merged==1))


#         case 'extensive'
# # extensive values dep on the size of the system
# # for example when downscaling the liquid water content of one
# # cell into two equally sized cells then the lwc in the new
# # cells are half of the lwc of the original ones
# # example: depth_old = 1:4 lwc_old = [1 0 0 1]
# #           depth_new = [0.1 0.2 0.6 1.2 3.5]

#             if depth_new()>depth_old()
#                 thick_last_old = depth_old()-depth_old(-1)
#                 depth_old = [depth_old, depth_new()]
#                 thick_last_old_new = depth_old()-depth_old(-1)
#                 var_old = [var_old, var_old()/thick_last_old*thick_last_old_new]


#             depth_merged = sort([depth_old depth_new])

#             thick_merged = depth_merged
#             thick_merged(2:) = depth_merged(2:) - depth_merged(1:-1)

#             thick_old = depth_old
#             thick_old(2:) = depth_old(2:) - depth_old(1:-1)

#             ind_bin = discretize(depth_merged,[0 depth_old],'IncludedEdge','right')

#             if sum(isnan(ind_bin))>1
#                 error('Some depths asked in depth_new not covered by depth_old')

#             var_merged = var_old(ind_bin).* thick_merged ./ thick_old(ind_bin)

#             ind_bin_new = discretize(depth_merged,[0 depth_new],'IncludedEdge','right')
#             var_merged(isnan(ind_bin_new)) =  []
#             ind_bin_new(isnan(ind_bin_new)) =  []

#             var_new = accumarray(ind_bin_new', var_merged')'


#     if transpose_at_the_==1
#        var_new=var_new'

# return var_new
