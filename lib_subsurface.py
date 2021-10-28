
def [ psnowc, psnic, pslwc, ptsoil, zrfrz, prhofirn,\
    zsupimp, pdgrain, zrogl, psn, pgrndc, pgrndd, pgrndcapc, pgrndhflx,\
    dH_comp, snowbkt, compaction, c] \
    = subsurface(pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, \
    ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep,\
    snowbkt, c)

# HIRHAM subsurface scheme - version 2016
# Developped by Peter Langen (DMI), Robert Fausto (GEUS)
# and Baptiste Vandecrux (DTU-GEUS)
#
# Variables are:
#     grndhflx - diffusive heat flux to top layer from below (W/m2, positive upwards)
#     tsoil - Layerwise temperatures, top is layer 1 (K)
#     slwc - Layerwise liquid water content (m weq)
#     snic - Layerwise ice content (m weq)
#     snowc - Layerwise snow content (m weq)
#     rhofirn - Layerwise density of snow (kg/m3)
#     slush - Amount of liquid in slush bucket (m weq)
#     snmel - Melting of snow and ice (daily sum, mm/day weq)
#     rogl - Runoff (daily sum, mm/day weq)
#     sn - Snow depth (m weq)
#     rfrz - Layerwise refreezing (daily sum, mm/day weq)
#     supimp - Superimposed ice formation (daily sum, mm/day weq)
#  
# thickness_act(n) = snowc(n)*(rho_w/rhofirn(n)) + snic*(rho_w/rho_ice) + slwc

# NEW Model set up:
# Each layer has a part which is snow (snowc), ice (snic) and water (slwc), 
# and the total water equivalent thickness of layer n is 
# thickness_weq(n) = snowc(n)+snic(n)+slwc(n)
# This thickness is allowed to vary within certain rules.

[ptsoil] = tsoil_diffusion (pts, pgrndc, pgrndd, ptsoil, c)


# pore_space = psnowc .* c.rho_water.*( 1./prhofirn - 1/c.rho_ice)
# excess_ice = max(0, psnic * c.rho_water / c.rho_ice - pore_space)
# thickness_act = psnowc.*(c.rho_water./prhofirn) + excess_ice
# depth = cumsum(thickness_act)
# # plot(T_firn_ini-273,-aux ,'-o')
# # hold on
# plot(ptsoil-273,-depth ,'-o')
# ylim([-5 0])
# xlim([-25 -15])
# hold off
# disp([pgrndc pgrndd])
# pause(0.1)
[prhofirn, dH_comp, compaction] = densification (pslwc, psnowc , prhofirn, ptsoil, c)
        if c.track_density
            [c.rhoCC20_aft_comp(1), c.rhoCC20_aft_comp(2)] = Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
        

[pdgrain] =  graingrowth ( pslwc, psnowc , pdgrain, c)

[psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]\
    = snowfall_new (zsn, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, pts,snowbkt, zraind, zsnmel,  c)
        # Update BV 2018: Density and temperature tracker
        if c.track_density
            [c.rhoCC20_aft_snow(1), c.rhoCC20_aft_snow(2)] = \
                Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
        

#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]\
    = sublimation_new (zsn, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, snowbkt, c)
        # Update BV 2018: Density and temperature tracker
        if c.track_density
            [c.rhoCC20_aft_subl(1), c.rhoCC20_aft_subl(2)]= Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
        
        
#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil ] \
    = rainfall_new (zraind, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, pts, c)

#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[zso_capa, zso_cond] = ice_heats (ptsoil, c)

[psnowc, psnic, pslwc, snowbkt] = \
        melting_new (psnowc, psnic, pslwc, zsnmel, snowbkt, ptsoil,prhofirn, c)

        # Update BV 2018: Density and temperature tracker
        if c.track_density
             [c.rhoCC20_aft_melt(1), c.rhoCC20_aft_melt(2)]= \
                 Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
                      

#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

if c.hetero_percol
    [pslwc] = hetero_percol (prhofirn, psnowc , psnic, pslwc, pdgrain, c)


[prhofirn, psnowc , psnic, pslwc , pdgrain, zrogl] =\
    perc_runoff_new (prhofirn, psnowc , psnic, pslwc , \
    pdgrain, c)

        # Update BV 2018: Density and temp tracker
        if c.track_density
            [c.rhoCC20_aft_runoff(1), c.rhoCC20_aft_runoff(2)] = \
                Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
        

#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[psnic, pslwc, ptsoil, zrfrz]\
    = refreeze (psnowc, psnic, pslwc, ptsoil, c )
        
#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[ptsoil ,  psnic, pslwc, zsupimp]\
    =  superimposedice (prhofirn, ptsoil            \
    , psnowc  , psnic, pslwc, zso_cond, c )

        # Update BV 2018: Density and temperature tracker
        if c.track_density
            [c.rhoCC20_aft_rfrz(1),c.rhoCC20_aft_rfrz(2)]  = Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
        

[prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain] =\
    merge_small_layers (prhofirn, psnowc , psnic, pslwc, ptsoil , \
    pdgrain, c)

        # Update BV 2018: Density and temperature tracker
        if c.track_density
            rho_avg_test = Calculate20mAvgDensity(psnowc, psnic, pslwc, snowbkt, prhofirn, ptsoil, c)
            if abs(rho_avg_test - c.rhoCC20_aft_rfrz(1)) > 1
                error('merge_small_layer changing density')
            
        

#BV 2017 updating cdel, cmid and rcdel
c = update_column_properties(c,psnowc, psnic, pslwc)

[psn] =  calc_snowdepth1D (psnowc, psnic, snowbkt, c)

[pgrndc, pgrndd, pgrndcapc, pgrndhflx]\
    = update_tempdiff_params (prhofirn, pTdeep \
    , psnowc, psnic, ptsoil, zso_cond, zso_capa, c)





def [prhofirn, dH_comp, dV] = densification (pslwc, psnowc , prhofirn, ptsoil, c)
#   densification: updates the subsurface firn/snow density with
#   gravity-driven densification.
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

#   ! PLA densification (Feb 2015)
#   ! Densification of layers (RSF+PLA Jan-Feb 2015):
#   ! Vionnet et al. (GMD, 2012) densification eqs 5-9

    #            ! Pressure in Pa (evaluated at layer mid-point).
    #            ! We have neglected surface slope in eqn 6.

    
    sigma_pres(1:c.jpgrnd,1) = c.cmid * c.rho_water * c.g 
    
    f1f2(1:c.jpgrnd,1) = 4./( 1 + 60.*(pslwc./ (psnowc + c.smallno)) \
        .* (prhofirn ./ c.rho_water) ) # ! f2 = 4
    
    eta_firn(1:c.jpgrnd,1)   = f1f2 .* 7.62237e+6 .* (c.a_dens * prhofirn / 250 ) \
        .* exp( 0.1.*(c.T_0 - ptsoil) + c.b_dens * 0.023.*prhofirn )
    prhofirn_old = prhofirn
    prhofirn = prhofirn + c.zdtime*prhofirn .* sigma_pres ./eta_firn

    # Update Baptiste 2018
#     MO = NaN(size(prhofirn))
#     MO(prhofirn <= 550)  = max(0.25, 1.042 - 0.0916 * log(c.accum_AWS))
#     MO(prhofirn > 550)  = max(0.25, 1.734 - 0.2039 * log(c.accum_AWS))
# 
#     C = NaN(size(prhofirn))
#     C(prhofirn <= 550)  = 0.07
#     C(prhofirn > 550)  = 0.03
#     
#     Ec = 60*1000 # Jmol?1
#     Eg = 42.4*1000 # J mol?1
#     prhofirn = prhofirn + c.zdtime.*\
#         MO .* C .* c.accum_AWS * c.rho_water/365/24/3600 .* c.g .* (c.rho_ice - prhofirn) \
#         .* exp(-Ec./c.R./ptsoil + Eg./c.R./ptsoil)

# Arthern
#     c.R = 8.314
#     Ec = 60*1000
#     Eg = 42.4*1000
#     m0 = 1 #1.042 - 0.0916*log(acc)
#     m1 = 1 #1.734 - 0.2039*log(acc)
#     # since the standard time unit in the code is second
#     # we need to convert accumulation rates to kg/s
#     c_0 = 0.07 * c.mean_accum*c.rho_water/365/24/3600*c.g .* exp( - Ec / c.R ./ptsoil + Eg / c.R / c.Tdeep)
#     c_1 = 0.03 * c.mean_accum*c.rho_water/365/24/3600*c.g .* exp( - Ec / c.R ./ptsoil + Eg / c.R / c.Tdeep)
#     
#     ind1 = ( prhofirn <= 550)
#     prhofirn(ind1) = prhofirn(ind1) + c.zdtime * m0.* c_0(ind1).*(c.rho_ice - prhofirn(ind1))
#     ind2 = prhofirn > 550
#     prhofirn(ind2) = prhofirn(ind2)  + c.zdtime * c_1(ind2) .* m1 .*(c.rho_ice - prhofirn(ind2))
    
    prhofirn(prhofirn > c.rho_ice)= c.rho_ice

    dV = psnowc * c.rho_water ./prhofirn_old - psnowc * c.rho_water ./prhofirn
    dH_comp = sum(dV) #in real m, here positive when compaction o


def [dgraindtF] = dgraindtF(dg,pl,ps, c)
#     ! The F in the name stands for "def"
#     ! Calculates d(grainsize-diameter) /dt (in mm/s) as in Hirashima (2010), using saturated
#     ! and unsaturated rates from Brun (1989), Tusima (1978) and Katsuhima (2009).

#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================
#     ! Mass liquid water content in #
L = pl./(ps+c.smallno).*100

#     ! Brun derives a def for grain growth that increases as L^3 (and divides by d^2).
#     ! Katushima uses this for L<=10#. Beyond 10 # Katushima uses the minimum value of the
#     ! Brun formuc.lation and the constant-divided-by-d^2 found by Tusima.
#     ! However, the Tusima-constant is lower than the L^3 def at 10# so that would lead
#     ! to a discontinuous drop when going beyond 10#.
#     ! Instead, we do as we also think Hirashima does "Therefore, to consider grain growth
#     ! in the water-saturated layer, the grain coarsening model of Tusima (1978) was used
#     ! as an upper boundary to saturated grain growth following Katsushima et al. (2009a)":

#     ! Take the L^3 def by Brun up until this becomes larger than the Tusima-constant:
Brun   = 2./pi.* ( 1.28e-8 + 4.22e-10 .* L.^3)
Tusima = 6.94e-8
#     ! d (grainsize-diamter) /dt in mm/s
dgraindtF = 1./(dg.^2).*min(Brun,Tusima)



def [pdgrain] =  graingrowth ( pslwc, psnowc , pdgrain, c)
#   graingrowth: Update the subsurface grain size vector after grain growth
#   due to metamorphism.

#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

#   ! PLA Darcy (may 2016)
pdgrain( 1:c.jpgrnd) = pdgrain( 1:c.jpgrnd) + \
    c.zdtime*dgraindtF(pdgrain( 1:c.jpgrnd),\
    pslwc( 1:c.jpgrnd),psnowc(1:c.jpgrnd), c)


def [zso_capa, zso_cond] = ice_heats (ptsoil, c)
#ice_heats:computes the subsurface thermal capacity [J/K] and thermal 
# conductivity zso_cond [J/S/M/K] from the subsurface temperature diffusivity
# c.zdifiz [M**2/S]
# Input:
#   ptsoil - subsufrace temperature in Kelvin
#   c - structure containing the constants
# Output:
#   zso_capa - Volumetric heat capacity of ice for each of the layers in
#   J/kg/K. Note that it is the same for snow!
#   zso_cond - Thermal conductivity of ice for each of the layer, in W/m/K
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

# Update BV2017
# Volumetric heat capacity of ice. In J/m^3/K.
zso_capa = zsn_capaF(c.rho_ice, ptsoil)

# update BV2017
# thermal conductivity of pure ice in W/m/K according to Yen 1981
zso_cond = 9.828*exp(-0.0057*ptsoil)

# old version
# zso_cond = zso_capa*c.zdifiz


def [psnowc, psnic, pslwc, snowbkt] = \
    melting_new (psnowc, psnic, pslwc, zsnmel, snowbkt, ptsoil, prhofirn, c)
# melting: Transform an amount zsnmel (m eq) of frozen material into liquid
# water. Starts at the surface and continues downward as long as more
# material needs to be melted. Update BV2017: The proportion of ice and snow  
# that is melted is given by the proportion of ice and snow present in the 
# layer. Material from the snowbkt is taken before anything else. The
# volume of what is being melted is calculated.
#	
#   input:
#         zsnmel - mass (m weq) of material that needs to be melted
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         snowbkt - Fresh snow bucket where new snow is added. When its
#         content is greater than lim_new_lay, then a new layer is created.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]
#           vol_melted - volume of the material that was melted
#
#   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
#   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
#=========================================================================

# First we handle melting. Make sure all melt energy is used. Start by melting in
# top layer. If necessary go to deeper layers.
zsnout = zsnmel

jk = 0

# first  we melt the content of the snow bucket
if (snowbkt > c.smallno && zsnout > c.smallno)
    zdel = min(snowbkt, zsnout)
    snowbkt = snowbkt - zdel
    zsnout = zsnout - zdel
    pslwc(1) = pslwc(1) + zdel


#         ! Enter if there is any melting
while (zsnout > c.smallno)
    #  Start at top and work downward if necessary
    jk = jk + 1 
    
    # Update BV2017
    # How much energy is needed to bring the layer to melting point
    deltaT = c.T_0 - ptsoil(jk)
    # volumetric heat capacity of the layer in J/m3/K
    heat_capa_vol = zsn_capaF(prhofirn(jk),ptsoil(jk))
    volume = (psnic(jk) + psnowc(jk))*c.rho_water / prhofirn(jk) 
    # heat capacity in J/K
    heat_capa = heat_capa_vol .* volume
    # energy needed to warm layer to 0degC in J
    warming_energy = deltaT*heat_capa
    if zsnout - warming_energy/c.L_fus/c.rho_water>0
        # if more melt than cold content
        # converting that energy into meter of meltwater and removing it from the melt
        zsnout = zsnout - warming_energy/c.L_fus/c.rho_water
    else
        # if less melt than cold content then melt energy is used to warm
        # the layer
        deltaT = zsnout/heat_capa*c.L_fus*c.rho_water
        ptsoil(jk) = ptsoil(jk) + deltaT
    

    #  How much frozen mass do we have in layer available for melting?       
    zdel = min(zsnout, psnowc(jk) + psnic(jk))

    # Update BV2017: Snow and ice are now melted simultaneously
    snow_melt = psnowc(jk)./(psnowc(jk)+psnic(jk)) .* zdel
    ice_melt = psnic(jk)./(psnowc(jk)+psnic(jk)) .* zdel
    
    psnowc(jk) = psnowc(jk) - snow_melt
    psnic(jk) = psnic(jk) - ice_melt
    pslwc(jk) = pslwc(jk) + zdel
    
    zsnout = zsnout - zdel


if sum((psnowc+psnic)==0)>1
    disp('MELTING MORE THAN ONE LAYER')
    disp(zsnmel)




def [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = \
            merge_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c)
# merge_layer: This def finds the two layers in the column that are
# most similar,merges them and frees the top layer. The similarity between
# layers is a weighted average of 5 objective criteria: difference in
# temperature, in firn denisty, in grain size, in ice fraction, in water
# content plus a 6th crietria that encourage the model to merge layers that
# are deeper in the column and therefore preserve good resolution near the
# surface. The parameter w can be tuned to encourage deep merging strongly(
# w<<1) or not encouraging at all (w>>1).
#	
#   input:
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil]
#
#   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
#   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
#=========================================================================

# CHANGE TO PERMEABILITY CRITERIA?
    depth_weq = psnowc + psnic
    delta_depth = depth_weq(2:)- depth_weq(1:-1)
    # 1st criterion: Difference in temperature
    diff = abs(ptsoil(1:-1) - ptsoil(2:))./delta_depth
    diff(+1) = diff()
    crit_1 = interp1([0, c.merge_crit_temp/0.1], [1, 0], diff, 'linear', 0)

    # 2nd criterion: Difference in firn density
    diff = abs(prhofirn(1:-1) - prhofirn(2:))./delta_depth
    diff(+1) = diff()
    crit_2 = interp1([0, c.merge_crit_dens/0.1], [1, 0], diff, 'linear', 0)
    crit_2(and(psnowc(1:-1)==0 , psnowc(2:)==0)) = 0

    # 3rd criterion: Difference in grain size
    diff = abs(pdgrain(1:-1) - pdgrain(2:))./delta_depth
    diff(+1) = diff()
    crit_3 = interp1([0, c.merge_crit_gsize/0.1], [1, 0], diff, 'linear', 0)
    crit_3(and(psnowc(1:-1)==0 , psnowc(2:)==0)) = 0

    # 4th criterion: Difference in water content (rel. to layer size)
    # saturation instead?
    rel_lwc = pslwc./(pslwc+psnic+psnowc)
    diff = abs(rel_lwc(1:-1) - rel_lwc(2:))./delta_depth
# 	diff(+1) = diff()
    crit_4 = interp1([0, c.merge_crit_rel_lwc/0.1], [1, 0], diff, 'linear', 0)
    crit_4(and(pslwc(1:-1)==0 , pslwc(2:)==0)) = 0
    crit_4(+1) = crit_4()

    # 5th criterion: difference in ice content (rel. to layer size)
    rel_ic = psnic./(pslwc+psnic+psnowc)
    diff = abs(rel_ic(1:-1) - rel_ic(2:))./delta_depth
# 	diff(+1) = diff()
    crit_5 = interp1([0, c.merge_crit_rel_snic], [1, 0], diff, 'linear', 0)
    crit_5(and(psnic(1:-1)==0 , psnic(2:)==0)) = 0
    crit_5(+1) = crit_5()

    # 6th criterion: arbitrary preference for deep layers to merge more
    # than shallow layers
    thickness_weq = psnic + psnowc + pslwc
    depth_weq = cumsum(thickness_weq)
    mid_point_weq = depth_weq - thickness_weq./2

    isshallow = mid_point_weq < c.shallow_firn_lim
    isdeep    = mid_point_weq > c.deep_firn_lim

    crit_6(isshallow) = 0
    crit_6(isdeep) = 1
    
    if sum(isshallow+isdeep) < length(mid_point_weq)
       isbetween = find(and(~isshallow,~isdeep))
        crit_6(isbetween) = interp1(\
            [c.deep_firn_lim, c.shallow_firn_lim],\
            [1, 0], mid_point_weq(isbetween))
    
    crit_6 = crit_6'

    # Update BV2020: max layer thickness
    crit_7=interp1([0 c.max_lay_thick],[1 -1],[thickness_weq(2:) thickness_weq()], 'linear', -1)
    
    # final rating:
    w1 = 1
    w2 = 1
    w3 = 1
    w4 = 1
    w5 = 1
    w6 = 2
    w7 = 3
    crit = (w1.*crit_1 + w2.*crit_2+ w3.*crit_3 +\
        w4.*crit_4 + w5.*crit_5 + w6 .* crit_6+ w7 .* crit_7)./(w1+w2+w3+w4+w5+w6+w7)
    i_merg = find(crit == max(crit),1,'last')
    i_merg = min(c.jpgrnd-1, i_merg)

# fprintf('Merging layer #i and #i.\n',i_merg, i_merg +1)
# subplot(2,1,1)
# hold off
# plot(crit_1)
# hold on
# plot(crit_2)
# plot(crit_3)
# plot(crit_4)
# plot(crit_5)
# plot(crit_6)
# plot(crit_7)
# plot(crit,'LineWidth',2)
# leg(num2str([1:8]'),'location','eastoutside')
# plot(i_merg*[1 1], [-1 1],'LineWidth',2)
# 
# subplot(2,1,2)
# hold off
# stairs0(depth_weq,prhofirn,'LineWidth',2)
# hold on
# for i = 1:length(depth_weq)
# stairs0(depth_weq(i)*[1 1],[0 800])
# 
# 
# pause(0.1)
    # layer of index i_merg and i_merg+1 are merged
    if (psnowc(i_merg+1) + psnowc(i_merg))> c.smallno
        if psnowc(i_merg+1)<c.smallno
            psnowc(i_merg+1) = 0
        
        if psnowc(i_merg)<c.smallno
            psnowc(i_merg) = 0
        
        #if there is snow in the two layers
        pdgrain(i_merg+1) = (psnowc(i_merg) * pdgrain(i_merg) + psnowc(i_merg+1) * pdgrain(i_merg+1)) \
           /(psnowc(i_merg+1) + psnowc(i_merg))
        prhofirn(i_merg+1) = (psnowc(i_merg+1) + psnowc(i_merg)) \
           / (psnowc(i_merg)/prhofirn(i_merg) +  psnowc(i_merg+1)/prhofirn(i_merg+1))
    else
        pdgrain(i_merg+1) = -999
        prhofirn(i_merg+1) = -999
    
    ptsoil(i_merg+1) = ((psnic(i_merg) + psnowc(i_merg) + pslwc(i_merg)) * ptsoil(i_merg) \
        + (psnic(i_merg+1) + psnowc(i_merg+1) + pslwc(i_merg+1)) * ptsoil(i_merg+1)) \
        /(psnic(i_merg) + psnowc(i_merg) + pslwc(i_merg) + psnic(i_merg+1) + psnowc(i_merg+1) + pslwc(i_merg+1))

    psnowc(i_merg+1) = psnowc(i_merg+1) + psnowc(i_merg)
    psnic(i_merg+1) = psnic(i_merg+1) + psnic(i_merg)
    pslwc(i_merg+1) = pslwc(i_merg+1) + pslwc(i_merg)

    # shifting column to free top layer
    if i_merg ~= 1
        psnowc(2:i_merg) = psnowc(1:i_merg-1)
        psnic(2:i_merg) = psnic(1:i_merg-1)
        pslwc(2:i_merg) = pslwc(1:i_merg-1)
        pdgrain(2:i_merg) = pdgrain(1:i_merg-1)
        prhofirn(2:i_merg) = prhofirn(1:i_merg-1)
        ptsoil(2:i_merg) = ptsoil(1:i_merg-1)
    


def [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain] =\
    merge_small_layers (prhofirn, psnowc , psnic, pslwc, ptsoil , \
        pdgrain, c)

# Now that the water has percolated and refrozen, we go through the column and check
# that there is no layer that ed up too small

for jk =1:c.jpgrnd-1

    too_small = (psnic(jk) + psnowc(jk) + pslwc(jk)) < c.lim_new_lay
    while too_small== 1
        if (psnic(jk) + psnowc(jk) + pslwc(jk)) > c.smallno
            # then we merge this layer with the layer below
            if (psnowc(jk+1) + psnowc(jk))~=0
                pdgrain(jk+1) = (psnowc(jk) * pdgrain(jk) + psnowc(jk+1) * pdgrain(jk+1)) \
                   /(psnowc(jk+1) + psnowc(jk))
                prhofirn(jk+1) = (psnowc(jk+1) + psnowc(jk)) \
                   / (psnowc(jk)/prhofirn(jk) +  psnowc(jk+1)/prhofirn(jk+1))
            else
                pdgrain(jk+1) = c.dgrainNew
                prhofirn(jk+1) = c.rho_ice
            
            ptsoil(jk+1) = ((psnic(jk) + psnowc(jk) + pslwc(jk)) * ptsoil(jk) \
                + (psnic(jk+1) + psnowc(jk+1) + pslwc(jk+1)) * ptsoil(jk+1)) \
                /(psnic(jk) + psnowc(jk) + pslwc(jk) + psnic(jk+1) + psnowc(jk+1) + pslwc(jk+1))

            psnowc(jk+1) = psnowc(jk+1) + psnowc(jk)
            psnic(jk+1) = psnic(jk+1) + psnic(jk)
            pslwc(jk+1) = pslwc(jk+1) + pslwc(jk)

            # now layer jk is free but we need to shift the column from 1 to jk
            # in order to free the first layer instead (just for being
            # compatible with the SplitLayer def)
            if jk>1
                psnowc(2:jk) = psnowc(1:jk-1)
                psnic(2:jk) = psnic(1:jk-1)
                pslwc(2:jk) = pslwc(1:jk-1)
                pdgrain(2:jk) = pdgrain(1:jk-1)
                prhofirn(2:jk) = prhofirn(1:jk-1)
                ptsoil(2:jk) = ptsoil(1:jk-1)
                # now top layer is free and another layer can be splitted
                # elsewhere
            
        

        [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = \
            split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c)

        # if a layer was split above indice jk then the new, merged layer
        # is now located at jk+1 and will be the next one tested in the
        # for loop. However if the split happened below indice jk then
        # the merged layer is still located at indice jk and ti should be
        # tested again to see if this merging was enough
        too_small = (psnic(jk) + psnowc(jk) + pslwc(jk)) < c.lim_new_lay
        
        # in the first case, the new 'too_small' will be the one of the layer
        # just above the one that was just merged and should be 0 because
        # the loop already went through that layer. In the second case, the
        # new 'too_small' is the one of the merged layer and if still too
        # small, would be set to 1 and another layer would be merged to it.
    



def [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil ] \
    = rainfall_new (zraind, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, pts, c)

#add_rainfall: Routine that adds rain water input on top of the 
#subsurface column. Update BV2017: No mass shift anymore. Water is just
#added to the water fraction of the top layer and the total mass of the layer is
#allowed to change.
#	
#   input:
#         zraind - amount (m weq) of input rain.
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         pts - surface temperature (K)
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil]
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

if zraind > c.smallno
    T_rain_in = max(pts,c.T_0)
    ptsoil(1) = ((psnowc(1) + psnic(1) + pslwc(1)) * ptsoil(1)  + zraind * T_rain_in) \
        / (psnowc(1) + psnic(1) + pslwc(1) + zraind)

    pslwc(1) = pslwc(1) + zraind
    # rain switch the snowbkt into the first layer to mimic the fresh snow
    # getting wet



def [psnic, pslwc, ptsoil, zrfrz]\
    = refreeze (psnowc, psnic, pslwc, ptsoil, c )
#refreeze: Calculates the amount of refreezing undergone in the subsurface
#column after the meltwater has percolated.
#   syntax:
#       [psnowc, psnic, pslwc, ptsoil, zrfrz, pts ]\
#           = refreeze (psnowc, psnic, pslwc, ptsoil, c)
#   input:
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (mm weq). Their sum
#         for each layer should always be equal to the layer fixed water
#         eq. thickness.
#
#         ptsoil - vector of subsurface temperature (K)
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated prhofirn, ptsoil, psnowc, psnic, pslwc
#          zrfrz - amount of water (mm eq) refrozen during the time step
#          pts - surface temperature
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to IDL by
#   Robert S. Fausto and later to Matlab by Baptiste Vandecrux
#   (bava@byg.dtu.dk).
# ========================================================================

#======================================================
#  Here we do the refreezing based on the cold content
# of each layer converting mass from liquid to ice.
#======================================================
zrfrz=zeros(c.jpgrnd,1)
zpotref=zeros(c.jpgrnd,1)
coldcontent=zeros(c.jpgrnd,1)
cfrozen=zeros(c.jpgrnd,1)

# Determine layer cold content and convert mass from liq to ice
cfrozen(1:c.jpgrnd) = psnowc(1:c.jpgrnd) + psnic(1:c.jpgrnd)
coldcontent(1:c.jpgrnd) = max(0, (c.T_0-ptsoil(1:c.jpgrnd)))

# update summer 2013 by BV
# only a fraction C_k * dt of the potential refrozen water is allowed to
# refreeze during the time step dt. See Illangaskare et al. 1990, eq. 13
zpotref(1:c.jpgrnd) = c.Ck * \
    coldcontent .* cpiceF(ptsoil) .* cfrozen / c.L_fus
zrfrz(1:c.jpgrnd)= min(zpotref , pslwc)  
    
pslwc(1:c.jpgrnd) =  pslwc(1:c.jpgrnd) - zrfrz(1:c.jpgrnd)
psnic = psnic + zrfrz
ptsoil_old = ptsoil
# Update temperature with c.latent heat of freezing
ptsoil(1:c.jpgrnd) = ptsoil(1:c.jpgrnd) + \
    zrfrz(1:c.jpgrnd).* c.L_fus./( cpiceF(ptsoil).*(psnic(1:c.jpgrnd)-zrfrz(1:c.jpgrnd) +\
    psnowc(1:c.jpgrnd)) )
# pts = ptsoil(1)


def [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]\
    = snowfall_new (zsn, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, pts, snowbkt, zraind, zsnmel, c)
# snowfall_new: Routine that adds new fallen snow (net from 
# sublimation)on top of the subsurface column. It is first accumulated in 
# snowbkt and only when snowbkt exceed a threshold then a new layer is
# created. Update BV2017: No mass shift anymore. When a new layer is created, two 
# layers at depth are merged.
#	
#   input:
#         zsn - amount (m weq) of new fallen snow net of sublimationaning
#         that it can be negative if sublimation is greater than
#         snow accumulation).
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         pts - surface temperature (K)
#
#         snowbkt - Fresh snow bucket where new snow is added. When its
#         content is greater than lim_new_lay, then a new layer is created.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]
#
#   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
#   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
#=========================================================================

if (zsn > 0) # ! zsn means snowfall in timestep

    #snowfall added to the fresh snow bucket
    snowbkt = snowbkt + zsn
    
    if snowbkt > c.lim_new_lay
        #enough material for a new layer           
#         if sum(psnic+psnowc+pslwc)> c.lim_max_column
#         # Update BV2017
#         # if the column (excluding last layer) is larger than a threshold, 
#         # then last layer is discarded, all the layers are shifted downward
#         # and the new snow layer is put in first layer.
#             psnowc(2:) = psnowc(1:-1)
#             psnic(2:) = psnic(1:-1)
#             pslwc(2:) = pslwc(1:-1)
#             pdgrain(2:) = pdgrain(1:-1)
#             prhofirn(2:) = prhofirn(1:-1)
#             ptsoil(2:) = ptsoil(1:-1)
#         else
            [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = \
                merge_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil,c)
#         
#         disp('SNOWFALL')
#         disp(snowbkt)
        psnowc(1) = snowbkt
        psnic(1) = 0
        pslwc(1) = 0

        ptsoil(1) = min(pts,c.T_0)
        prhofirn(1) = c.rho_fresh_snow
        pdgrain(1) = c.dgrainNew
        snowbkt = 0
    


# Update BV 2017
# if there is rain or melt, the fresh snow bucket is added to the snow
# compartment of the first layer
if and(or(zraind>c.smallno,zsnmel>c.smallno),snowbkt>c.smallno)
    # mass-weighted average for temperature
    ptsoil(1) = (ptsoil(1)*(psnowc(1)+pslwc(1)+psnic(1)) + min(pts,c.T_0)*snowbkt)\
        /(psnowc(1)+pslwc(1)+psnic(1)+ snowbkt)

    #snow-mass-weighted average for grain size
    pdgrain(1) = (psnowc(1) * pdgrain(1) + snowbkt * c.dgrainNew) \
           /(psnowc(1) + snowbkt)
    
    # volume-weighted average for density
    prhofirn(1) = (psnowc(1) + snowbkt) \
       / (snowbkt/c.rho_fresh_snow +  psnowc(1)/prhofirn(1))
 
   psnowc(1) = psnowc(1) + snowbkt
   snowbkt = 0       

    


def [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = \
            split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c)
# split_layer: This def finds the layer that is the best to be
# splitted. In most cases, since gradients in temperature, water content
# and density are greater close to the surface, it will select the
# uppermost layer that has enough mass to be splitted (mass >
# lim_new_lay*2). If all layers in the column are too small to be split,
# then the model will create a new layer at the bottom of the column with
# 10 m of ice and a user defined temperature of Tdeep.
# WARNING: the top layer needs to be free: its content will be discarded
# during the renumbering of the column.
#	
#   input:
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil]
#
#   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
#   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
#=========================================================================   

    thickness_weq = psnic + psnowc + pslwc
    
    if and(sum(thickness_weq > c.lim_new_lay*2) > 0,\
            sum(thickness_weq) > c.min_tot_thick)
#         if there are some layers that have are thick enough to be
#         splitted into two layers of which thickness would be higher than
#         the limit thickness for new layer c.lim_new_lay

        depth_weq = cumsum(thickness_weq)
        mid_point_weq = depth_weq - thickness_weq./2

        #first criterion: Layer that are thicker are more likely to be splitted
        crit_1 = interp1([0 c.lim_new_lay  c.thick_crit], [0 0 1], thickness_weq,'linear', 1)

        #second criterion: Shallow layers are more likely to be splitted    
        crit_2 = interp1([0 c.deep_firn_lim], [1 0], mid_point_weq,'linear', 0)

        crit = (0.6*crit_1 + 0.4*crit_2) / 2
        
        #all layer that are too thin to be splitted are taken out
        too_thin = thickness_weq./2 < c.lim_new_lay
        crit(too_thin) = -1
        
        crit(1) = 0 #We don't want to split the first layer anyway
        [~, i_split] = max(crit)
#         fprintf('Splitting layer #i\n',i_split)

        psnowc(1:i_split-1) = psnowc(2:i_split)
        psnic(1:i_split-1) = psnic(2:i_split)
        pslwc(1:i_split-1) = pslwc(2:i_split)
        pdgrain(1:i_split-1) = pdgrain(2:i_split)
        prhofirn(1:i_split-1) = prhofirn(2:i_split)
        ptsoil(1:i_split-1) = ptsoil(2:i_split)
        # now layer i_split and i_split-1 are the same
        
        thickness_new = min(c.max_lay_thick, thickness_weq(i_split)/2)
        thickness_left = thickness_weq(i_split) - thickness_new
        r_new = thickness_new / thickness_weq(i_split)
        r_left = thickness_left / thickness_weq(i_split)

        psnowc(i_split-1) = r_new*psnowc(i_split)
        psnic(i_split-1) = r_new*psnic(i_split)
        pslwc(i_split-1) = r_new*pslwc(i_split)
        # prhofirn, ptsoil and pdgrain were already identical

        psnowc(i_split) = r_left*psnowc(i_split)
        psnic(i_split) = r_left*psnic(i_split)
        pslwc(i_split) = r_left*pslwc(i_split)
        # now layer is split
    else
#       if all layers are too small, then ice is taken from below.
        psnowc(1:-1) = psnowc(2:)
        psnic(1:-1) = psnic(2:)
        pslwc(1:-1) = pslwc(2:)
        pdgrain(1:-1) = pdgrain(2:)
        prhofirn(1:-1) = prhofirn(2:)
        ptsoil(1:-1) = ptsoil(2:)
        
        psnowc() = 0
        psnic() = c.new_bottom_lay
        pslwc() = 0
        pdgrain() = 1
        prhofirn() = 1
        ptsoil() = c.Tdeep
#         disp('ADDING NEW LAYER')
    



def [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]\
    = sublimation_new (zsn, psnowc, psnic, pslwc, pdgrain \
    , prhofirn, ptsoil, snowbkt, c)
# sublimation_new: Routine that removes sublimation from the first layer of
# the column. Update BV2017: The proportion of ice and snow that is 
# sublimated is given by the proportion of ice
# and snow present in the layer. If the first layer
# becomes smaller than lim_old_lay, it is merged with the second
# layer,another layer is splitted elsewhere and sublimation continues on 
# this new merged layer. Material from the snowbkt is taken before anything
# else.
#	
#   input:
#         zsn - mass (m weq) of material that needs to be sublimated
#         
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         pts - surface temperature (K)
#
#         snowbkt - Fresh snow bucket where new snow is added. When its
#         content is greater than lim_new_lay, then a new layer is created.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]
#
#   This script was developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
#   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
#=========================================================================

if (zsn < 0)  # zsn negative means net upward snowfall+sublimation
    zsnout = abs(zsn)  # zsnout is the mass the leaves column upwards
    
    # first  we sublimate the content of the snow bucket
    if snowbkt > c.smallno
        pot_subl = min(snowbkt, zsnout)
        snowbkt = snowbkt - pot_subl
        zsnout = zsnout - pot_subl
    
    
    #then comes the upper  layers of the column
    while (zsnout > c.smallno )

        # Update BV2017: Snow and ice are now melted simultaneously
        zdel = min(zsnout, psnowc(1) + psnic(1))
        snow_melt = psnowc(1)/(psnowc(1)+psnic(1)) * zdel
        ice_melt = psnic(1)/(psnowc(1)+psnic(1)) * zdel

        psnowc(1) = psnowc(1) - snow_melt
        psnic(1) = psnic(1) - ice_melt
        zsnout =zsnout - zdel

        if psnowc(1) + psnic(1) + pslwc(1) < c.lim_new_lay
            if (psnowc(2) + psnowc(1))> c.smallno
                #then merging layer 1 and 2 and placing it in layer 2
                pdgrain(2) = (psnowc(1) * pdgrain(1) + psnowc(2) * pdgrain(2)) \
                   /(psnowc(2) + psnowc(1))
                prhofirn(2) = (psnowc(2) + psnowc(1)) \
                   / (psnowc(1)/prhofirn(1) +  psnowc(2)/prhofirn(2))
#                Update PLA
                ptsoil(2) = ((psnic(1) + psnowc(1) + pslwc(1)) * ptsoil(1) \
                    + (psnic(2) + psnowc(2) + pslwc(2)) * ptsoil(2)) \
                    /(psnic(1) + psnowc(1) + pslwc(1) + psnic(2) + psnowc(2) + pslwc(2))


                psnowc(2) = psnowc(2) + psnowc(1)
                psnic(2) = psnic(2) + psnic(1)
                pslwc(2) = pslwc(2) + pslwc(1)
            else
                pdgrain(2) = -999
                prhofirn(2) = -999
            
            
            # now top layer is free and another layer can be splitted
            # elsewhere
            [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = \
                split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c)
        
    



def [cpice] = cpiceF(ptsoil)
#cpice: right now constant heat capacity of ice. Can be modified to be
#depant on ice caracteristics.
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================
# Specific heat capacity of ice J/kg/K
# cpice = 2.108e+03      # standard value at 0degC 

# Update BV2017
cpice = (2.7442 + 0.1282 .* ptsoil)/18*1000
# Yen 1981 originally in J/mol/K converted to J/kg/K
# 1kg of ice has same capacity as 1kg snow.



def [ptsoil , psnic, pslwc, zsupimp]\
    =  superimposedice (prhofirn, ptsoil            \
    , psnowc  , psnic, pslwc, zso_cond, c )
#superimposedice: Calculates the amount (zsupimp in mm eq) of liquid water 
# that is refrozen on top of perched ice layers. New version from 2016. 
# Does not use slush bucket.
#   syntax:
# [prhofirn, ptsoil , psnowc  , psnic, pslwc, zsupimp]\
#     =  superimposedice (prhofirn, ptsoil            \
#     , psnowc  , psnic, pslwc, zso_cond, c )
#	
#   input:
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (mm weq). Their sum
#         for each layer should always be equal to the layer fixed water
#         eq. thickness.
#
#         ptsoil - vector of subsurface temperature (K)
#
#         zso_cond - Subsurface thermal conductivity. See ice_heats
#         def.
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [psnowc , psnic, pslwc, ptsoil]
#          zsupimp - total amount of formed superimposed ice (mm weq)
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

# --------------------------------------
# ----- SUPERIMPOSED ICE FORMATION -----
# --------------------------------------

# In the "Darcy-version", superimposed ice formation is changed in these respects:
#   1 Water used is just slwc (and not a slush bucket)
#   2 This simplifies things because no mass is added, subtracted or shifted.
#     Just converted from slwc to snic.
#   3 The new ice used to be added to the cold ice layer below. Now it is simply
#     converted to ice in the layer, where the water resides (giving energy to
#     the cold layer below).
#   4 SIF can oc.ccur several places in the column, as there is no c.longer a slush
#     bucket. slwc is now allowed to be larger than potret (the excess just runs off
#     with the Zuo&Oerlemans timescale as the sluch water used to) and it is this
#     stock of slwc that is used for SIF.

# No SIF allowed onto inifinte sublayer, so only
# search for suitable water layer in 1,c.jpgrnd-1 (ice layer below is in 2,c.jpgrnd):
# jk is the layer our water is in. We then look at jk+1 to see if
# we can superimpose:
# PLA-BV 2016: volume (instead of mass) weighted average of density
rho_next(1:c.jpgrnd-1) = (psnic(2:c.jpgrnd)+psnowc(2:c.jpgrnd))./\
    ( psnic(2:c.jpgrnd)/c.rho_ice + psnowc(2:c.jpgrnd) ./ prhofirn(2:c.jpgrnd) )

nextisfrozen(1:c.jpgrnd-1) = (rho_next(1:c.jpgrnd-1) >= c.rho_pco )
isfrozen(1)=0
isfrozen(2:c.jpgrnd) = nextisfrozen(1:c.jpgrnd-1)
# The next layer has bulk density high enough that we will consider it as
# an ice layer. Therefore we calcuc.late SI-formation:

snowV = 0.5*psnowc(isfrozen==1)  .*c.rho_water./prhofirn(isfrozen==1)
iceV  = 0.5*psnic(isfrozen==1)  *c.rh2oice
totalV = snowV + iceV
zx1 = snowV ./ totalV
zx2  = iceV   ./ totalV
# update BV2017 zsn_cond is now vector
ki= 1./( zx1./zsn_condF(prhofirn(isfrozen==1), c) + zx2./zso_cond(isfrozen==1) )
dTdz =  (c.T_0-ptsoil(isfrozen==1))  ./ totalV 

# The potential SIformation (Only interested in positive SIF):
potSIform = max(0, ki.*dTdz / (c.rho_water*c.L_fus) * c.zdtime) 

# Make as much SI as there is water available in the layer above the ice (jk)
SIform = min(pslwc(nextisfrozen==1),potSIform)
zsupimp = zeros(size(pslwc))
zsupimp(nextisfrozen==1) = SIform

# Convert water to ice
pslwc(nextisfrozen==1) = pslwc(nextisfrozen==1) - SIform
#BV 2016 corrected to avoid underflow
psnic(nextisfrozen==1) = psnic(nextisfrozen==1) + SIform
 
# Update temperature of ice layer
ptsoil(isfrozen==1) = ptsoil(isfrozen==1) + \
    SIform * c.L_fus ./ (c.cdel(isfrozen==1).*cpiceF(ptsoil(isfrozen==1)))


def [ptsoil] = tsoil_diffusion (pts, pgrndc, pgrndd, ptsoil, c)
#   tsoil_diffusion: Update subsurface temperatures ptsoil based on previous 
#   time step coefficients for heat conduction
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

#         ! Upper layer
ptsoil(1)= pts

#         ! Lower layers
# # update BV2017
# for jk = 2:c.jpgrnd
#     ptsoil(jk) = pgrndc(jk) + pgrndd(jk) * ptsoil(jk)
# 
# WHY??
# original 
for jk = 1:c.jpgrnd-1
    ptsoil(jk+1) = pgrndc(jk) + pgrndd(jk) * ptsoil(jk)



def c = update_column_properties(c,psnowc, psnic, pslwc)

c.cdel=psnowc+pslwc+psnic
aux = cumsum(c.cdel)

c.cdelsum = 0

for jk = 1:c.jpgrnd
    c.cdelsum = c.cdelsum + c.cdel(jk)
    c.cmid(jk) = c.cdelsum - ( c.cdel(jk) / 2 )
    c.rcdel(jk) = 1/c.cdel(jk)




def [pgrndc,pgrndd,   pgrndcapc, pgrndhflx ]\
    = update_tempdiff_params (prhofirn, pTdeep                    \
    , psnowc, psnic, ptsoil, zso_cond, zso_capa, c)
# update_tempdiff_params: Updates the thermal capacity and conductivity of 
# subsurface column based on the new density and temperature profiles. Also
# calcualates subsurface heat flux to the surface pgrndhflx and calorific capacity of the
# ground pgrndcapc.
# Input:
#   zso_cond - Thermal conductivity of pure ice. Calculated in ice_heats.
#   zso_capa - Specific heat capacity of ice and snow.
#   pgrndc - heat capacity (W/K)
# 
# This script was originally developped by Peter Langen (pla@dmi.dk) and
# Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
# Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================


# PETER AND RUTH's VERSION (Note: we ignore the liquid content in all of 
# the following). We need physical layer thicknesses (i.e., in
# snow and ice meters, not liquid meters):

# We consider that each of our temperature in the column are taken at
# the center of the cell. In other words they are seen as nodes and not 
# as cells anymore. We then calculate the heat capacity and thermal
# conductivity of the material that separates neighbouring nodes.

snowV1=zeros(c.jpgrnd,1)
snowV2=zeros(c.jpgrnd,1)
iceV=zeros(c.jpgrnd,1)
pgrndc=zeros(c.jpgrnd,1)
pgrndd=zeros(c.jpgrnd,1)


## Calculate midpoint volume-weighted versions of capa (heat capacity) and
# kappa (thermal conductivity)

# PLA densification (Feb 2015) Update BV2017
# The first layer should be in equilibrium with the surface temperature.
# Therefore first two nodes are seperated by all of layer 1 and half of
# layer 2.
# volume of snow in upper layer:
snowV1(1) =  psnowc(1)*c.rho_water/prhofirn(1) 
# volume of snow in lower half layer:
snowV2(1) = 0.5*psnowc(2)*c.rho_water/prhofirn(2) 
# volume of ice in both layer:
iceV(1)  = (psnic(1) + 0.5*psnic(2)) *c.rho_water / c.rho_ice 

# For following nodes, two neighboring nodes are separated by half of the 
# upper layer and half of the lower layer
# volume of snow in upper half layer:
snowV1(2:c.jpgrnd) = 0.5*psnowc(2:c.jpgrnd)  .*\
    c.rho_water./prhofirn(2:c.jpgrnd)
# volume of snow in lower half layer:
snowV2(2:c.jpgrnd - 1) = 0.5*psnowc(3:c.jpgrnd)*c.rho_water./prhofirn(3:c.jpgrnd)
# volume of ice in both half layers:
iceV(2:c.jpgrnd - 1)  = 0.5*(psnic(2:c.jpgrnd - 1) + psnic(3:c.jpgrnd ))\
    .*c.rh2oice

# Bottom layer zcapa assuming below is ice to same thickness (at least)
snowV2(c.jpgrnd) = 0
iceV(c.jpgrnd) = 0.5 * (psnic(c.jpgrnd) + c.cdel(c.jpgrnd)) *c.rh2oice

# total mass separating two nodes
totalV = snowV1 + snowV2 + iceV

# ice and snow volumetric fractions in the material separating two nodes
snow_frac_lay_1 = snowV1./ totalV
snow_frac_lay_2 = snowV2./ totalV
ice_frac = iceV ./ totalV

# heat capacity of the layer calculated as the volume-weighted average of
# the snow and ice heat capacity. Mind the repeated index for the last
# layer. Here zcapa is still volumetric since everything on the right hand
# side has been deivided by 'totalV'. It is in J/m^3/K.
zcapa = snow_frac_lay_1 .* zsn_capaF(prhofirn, ptsoil) +\
    snow_frac_lay_2 .* \
    zsn_capaF(prhofirn([2:c.jpgrnd, c.jpgrnd]), ptsoil([2:c.jpgrnd, c.jpgrnd])) \
    + ice_frac .* zso_capa

# thermal conductivity of the layer calculated as the inverse of the sum of
# the inversed conductivities.
# thick_tt/k_tt_eff = thick_1/k_1 + thick_2/k_2 + thick_3/k_3
# Warning: It dos not include thermal exchange through air in pore and
# through water. In W/m/K.
zkappa = 1./ \
    ( snow_frac_lay_1./zsn_condF(prhofirn, c) + \
    snow_frac_lay_2./zsn_condF(prhofirn([2:c.jpgrnd, c.jpgrnd]), c) + \
    ice_frac./zso_cond )

## Calculate volumetric heat capacity and effective thermal conductivity
# by multiplying, resp. dividing, by each layer thickness to get to the 
# effective values

# calculating layer volume (without liquid water)
c.cdelV(1:c.jpgrnd) = psnowc(1:c.jpgrnd) *c .rho_water./prhofirn(1:c.jpgrnd) \
    + psnic(1:c.jpgrnd)*c.rh2oice
# !!! Why here not totalV ?

# Total heat capacity for each layer in W/K
zcapa_abs =zcapa .*c.cdelV ./ c.zdtime

# The following used to go jk=1,c.jpgrnd-1
# Now we include c.jpgrnd, because it is needed below:

# Real distance between midpoints. Note repeated index for last layer.
cmidV =  (c.cdelV([2:c.jpgrnd, c.jpgrnd])+c.cdelV)./2
# calculating effective thermal conductivity
zkappa_abs = zkappa ./ cmidV
# !!! why her not *cmidV

# We have now calculated our new versions of zdz1, zdz2, zcapa and zkappa.
# Before introducing diffusion with sub-model layer, the old code was used 
# as it were. Now, we for some more:

# In the original version, this loop calculated c and d of c.jpgrnd-1 
# (which are in turn used above (in next time step) to prognose 
# tsoil(c.jpgrnd). Now we calculate c and d of c.jpgrnd.
# These are not used to prognose t of c.jpgrnd+1 (because this isn't 
# relevant) but to initialize the upwards calculation of c and d:

# Sublayer is all ice (zcapa = zso_capa) and of mass c.cdel(c.jpgrnd)
# giving physical thickness c.cdel(c.jpgrnd)*c.rh2oice:
# Update BV2017: using the temperature depant heat capacity of ice
zcapa_abs_sublayer = zsn_capaF(c.rho_ice,pTdeep) * c.cdel(c.jpgrnd) \
    * c.rho_water / c.rho_ice / c.zdtime
# corresponds to zcapa_abs(c.jpgrnd+1)

z1 = zcapa_abs_sublayer + zkappa_abs(c.jpgrnd)

pgrndc(c.jpgrnd) = zcapa_abs_sublayer .* pTdeep ./ z1
pgrndd(c.jpgrnd)=zkappa_abs(c.jpgrnd)./z1
# PLA Tdeep (feb 2015) pTdeep changed in the above

# This loop went jk=c.jpgrnd-1,2,-1 (ie, calculating c and d for 3,2,1)
# Now, it goes   jk=c.jpgrnd,2,-1 (ie, calculating c and d for 4,3,2,1)
# It thus needs zdz1(c.jpgrnd) which is now also calculated above

for jk=c.jpgrnd:-1:2
    z1= 1/(zcapa_abs(jk) + zkappa_abs(jk-1) +  zkappa_abs(jk)*(1-pgrndd(jk)))

    pgrndc(jk-1)=(ptsoil(jk)*zcapa_abs(jk) +  zkappa_abs(jk)*pgrndc(jk))*z1
    
    pgrndd(jk-1)=zkappa_abs(jk-1)*z1

#   ---------------------------------------
#   COMPUTATION OFSIVE FLUX FROM GROUND AND
#   CALORIFIC CAPACITY OF THE GROUND:
#   ---------------------------------------------------------

pgrndhflx=zkappa_abs(1) * (pgrndc(1)          \
    +(pgrndd(1)-1) * ptsoil(1))
pgrndcapc=(zcapa_abs(1) * c.zdtime+                \
    c.zdtime * (1-pgrndd(1)) * zkappa_abs(1))


def [zsn_capaF] = zsn_capaF(RHOS,ptsoil)
# The F in the name zsn_capaF stands for "def"
zsn_capaF = cpiceF(ptsoil).*RHOS  # snow vol. heat capacity   [J/m**3/K]
# Here the specific (per kg) heat is converted into volumetric heat (per
# m^3) capacity by multiplying by the density.



# PLA densification (Feb 2015
def [zsn_condF] = zsn_condF(prhofirn, c)
# The F in the name zsn_condF stands for "def"

# update BV2017
zsn_condF = 2.22362*(prhofirn./1000).^1.88
# zsn_condF2 = cpice(0)*c.rho_ice*c.zdifiz*((prhofirn/c.rho_water).^1.88)
zsn_condF2 = (1-(c.rho_ice - prhofirn)/c.rho_ice)*2.1

# snow thermal conductivity [J/s/m/K]
# Yen, Y. (1981), Review of thermal properties of snow, ice and sea ice,
# Rep. 81-10, U. S. Army Cold Reg. Res. and Eng. Lab. (CRREL), Hanover, N. H


def [psn] =  calc_snowdepth1D (psnowc, psnic, snowbkt, c)
# calc_snowdepth1D: Diagnose snow depth (for use in albefor 
# parameterizations and other places?). Include only the snow above first
# perched ice layer (and include the snow of this layer)
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================
notice = 1 # Should be read as "this layer is not ice"

# Update BV2017: including snowbkt
if ( snowbkt < c.smallno )
    psn = 0
    # here we allow to look at the rest of the layers
else
    psn = snowbkt


jk = 1
# Then the next layers
while and(notice == 1, jk <= c.jpgrnd)
    if (psnic(jk) > c.icemax*c.cdel)
        #then it is considered as ice layer
        notice = 0
    else
        psn = psn + psnowc(jk)
    
    jk = jk+1




