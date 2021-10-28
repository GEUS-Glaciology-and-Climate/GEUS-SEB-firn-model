def [prhofirn, psnowc , psnic, pslwc, pdgrain, zrogl] =\
    perc_runoff_new (prhofirn, psnowc , psnic, pslwc, \
        pdgrain, c)

# perc_runoff_new: Calculates meltwater percolation and runoff in the column
# either according to a standard bucket scheme or to a Darcy-like bucket
# scheme. Update BV2017: no mass shift anymore, the water is
# moving from one layer to the next and layer's total mass is allowed to change.
# At the  of the routine, it is checked that no layer has a mass lower
# than lim_old_lay. If one is found too small, it is merged with its
# underlying layer and another layer is split elsewhere.
#
#   syntax:
#   [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain, pTdeep, zrogl] =\
#     perc_runoff_new (prhofirn, psnowc , psnic, pslwc, ptsoil , \
#           pdgrain, pTdeep, c.ElevGrad, zrogl, c)
#	
#   input:
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (m weq).
#
#         ptsoil - vector of subsurface temperature (K)
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         pTdeep - lower boundary temperature (K). Defined by the constant
#         T_ice_AWS so far.
#
#         zrogl - total amount of run off (mm weq)
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain, zrogl]
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#=========================================================================

#======================================================
#Here we do liquid water percolation and runoff
#======================================================
# *** Calculate potential Darcy fluxes for all columns ***
#  (if do-no-darcy, then fluxes will be calcuc.lated on-the-fly in loops below)
# This calculation makes sure potential darcy flux across interfaces does
# not exceed the water amount in the upper of the two layers

if (~c.do_no_darcy)
    [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic, pdgrain \
        , prhofirn, c)


#Update BV2017: t_runoff is now calculated outside of the subsurface scheme

# Bottom layer: Receive from above. Give excess liquid to runoff.
jk = c.jpgrnd
zrogl = 0
    # Remove runoff from layer (and add to runoff-box zrogl)
    # Adjust bottom layer interface due to removal of mass
    # PLA Darcy 2016
     if ( c.calc_CLliq )
             liqmaxloc = CLliqF(prhofirn(jk),c)
          else
             liqmaxloc = c.liqmax
     

    liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1)
    potret    = max( liqmaxM* psnowc(jk) , 0 )
    # Calculate liqexcess from layer water. Make sure it is not negative:
    liqexcess = max ( pslwc(jk) - potret , 0 )
    liqro = liqexcess/c.t_runoff * c.zdtime
    # Take runoff from water content and give to runoff box (Update PLA)
    zrogl    = zrogl    + liqro  
    pslwc(jk) = pslwc(jk) - liqro

for jk = c.jpgrnd-1:-1:1
        # BV2017 removing percolation blocking criteria
        if ThetaF(pslwc(jk+1), psnowc(jk+1), prhofirn(jk+1), c) >= 1
            # if next layer already saturated, then we do not allow flux 
            # into that layer, but instead activate runoff
            dflux = 0
            do_runoff = 1 
        else
            if ( c.do_no_darcy )
                # Potential Darcy-like flux in case of do-no-Darcy,
                # basically just all liqexcess:
                # PLA Darcy 2016
                 if ( c.calc_CLliq )
                         liqmaxloc = CLliqF(prhofirn(jk),c)
                      else
                         liqmaxloc = c.liqmax
                 

                liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1)
                potret    = max( liqmaxM* psnowc(jk) , 0 )
                liqexcess = pslwc(jk) - potret
                darcy_fluxes(jk) = max(liqexcess , 0 )
            

            # Calculate water in next layer, when this is at saturation (Theta = 1):
            plsat = psnowc(jk+1)*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk+1) - 1)

            # Do not allow flux to be greater than plsat-pl in next layer.
            # Also, do not allow it to be negative
            dflux = max ( min( darcy_fluxes(jk) , plsat-pslwc(jk+1) ) , 0  )

                # Update BV2017: if no snow to hold liquid water then runoff is
                # turned on
            if or( darcy_fluxes(jk) >= plsat-pslwc(jk+1),\
                    psnowc(jk)<c.smallno)
                # There is enough darcy flow to fill up layer below (and perhaps more)
                # so runoff should be activated (on what's left here after
                # we have moved down enough to fill layer below):
                do_runoff = 1
            else
                # No runoff from this layer, since darcy flow is not large enough
                # to fill up layer below
                do_runoff = 0
            
        

        if ( dflux > 0 )
            # Yes: Darcy flow to next layer
               # Update BV2017: Now temperature is not updated. water
               # percolating is at 0degC cannot change temperature until phase
               # change occurs in refreezing scheme
            pslwc(jk+1) = pslwc(jk+1) + dflux
            pslwc(jk) = pslwc(jk) - dflux 
        
        
        if and(do_runoff, ~c.avoid_runoff)
            # Yes: Remove runoff from layer (and add to runoff-box zrogl)
            # PLA Darcy 2016
             if ( c.calc_CLliq )
                     liqmaxloc = CLliqF(prhofirn(jk),c)
                  else
                     liqmaxloc = c.liqmax
             

            liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1)
            potret    = max( liqmaxM* psnowc(jk) , 0)
            
            # Calculate liqexcess from layer water. Make sure it is not negative:
                # Update BV2017: since dflux has already left slwc(jk) so no
                # need to remove it anymore (compared to old version)
            liqexcess = max ( pslwc(jk) - potret , 0 )
                       
                #Update BV2017: for all layer, if there is no snow to hold
                #liquid water then runoff is instantaneous ( not allowing
                #subglacial lakes to form
            if (psnowc(jk) < c.smallno)
                # Here in layer: If there is no snow, run off immediately
                # in other word surface runoff is instantaneous
                liqro = liqexcess
            else
                Theta = ThetaF( pslwc(jk),psnowc(jk),prhofirn(jk), c)
#                 liqro_darcy = kF( Theta,pdgrain(jk),prhofirn(jk),psnic(jk) ,psnowc(jk), c) * c.ElevGrad 
                # old version based on Zuo and Oerlemans (1996)
                liqro_darcy  = liqexcess/c.t_runoff * c.zdtime
                
                liqro = min(liqro_darcy,liqexcess)
                pore_space = psnowc(jk)*c.rho_water*(1/prhofirn(jk)-1/c.rho_ice)

                #Update PLA
                if pslwc(jk) - liqro > pore_space
                    liqro = pslwc(jk) - pore_space
                
                    
#                 fprintf('Runoff #f #f #f\n',liqro_darcy, liqro_zo, liqro)
            
            
            # Take runoff from water content and give to runoff box
            zrogl    = zrogl    + liqro  
            pslwc(jk) = pslwc(jk) - liqro
        




def [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic,\
    pdgrain , prhofirn, c)

#calc_darcy_fluxes: Calculates the amount of water (mm weq) that each layer
#can transmit to the next one according to a Darcy flow.
#
#   syntax:
#   [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic,\
#            pdgrain , prhofirn, c)
#	
#   input:
#         psnowc, psnic, pslwc - vectors of respectively snow, ice and
#         liquid water part for each subsurface layer (mm weq). Their sum
#         for each layer should always be equal to the layer fixed water
#         eq. thickness.
#
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         prhofirn - vector of firn (snow) density (kg/m^3) 
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#         darcy_fluxes - vector containing the amount of water (mm weq) that each layer
#         can transmit to the next one according to a Darcy flow.
#
#   This script was originally developped by Peter Langen (pla@dmi.dk) and
#   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
#   Baptiste Vandecrux (bava@byg.dtu.dk).
#
#   Aug. 2016
#=========================================================================

darcy_fluxes = zeros(c.jpgrnd,1)
for jk = 1:c.jpgrnd-1
    # Special treatment of top layer in case of no snow
    if (jk == 1 && (psnowc(jk) < c.smallno))
        # In case of no snow, the calculation below of K crashes. Instead we say
        # * if layer below may receive water, it gets all there is in layer 1
        # * if layer may not receive water, it all runs off.
        #
        # As long as we set darcy_fluxes(1) = pslwc(1) then this will
        # be taken care of in the perc_runoff code.
        darcy_fluxes(jk) = pslwc(jk)
    else
        # Other layers do (should!) not come into a situation where there is only ice and water
        qlim = qlimF(                      \
            pslwc(jk)   , pslwc(jk+1)    , \
            psnowc(jk)  , psnowc(jk+1)   , \
            psnic(jk)   , psnic(jk+1)    , \
            prhofirn(jk), prhofirn(jk+1) , \
            pdgrain(jk) , pdgrain(jk+1), c )
        Theta1 = ThetaF( pslwc(jk)  ,psnowc(jk)   ,prhofirn(jk), c )
        Theta2 = ThetaF( pslwc(jk+1),psnowc(jk+1) ,prhofirn(jk+1), c)
        h1     = hHirF( Theta1, pdgrain(jk), c)
        h2     = hHirF( Theta2, pdgrain(jk+1), c)

        delz1 = c.rho_water*(psnowc(jk)  /prhofirn(jk)   + psnic(jk)  /c.rho_ice) # layer 1
        delz2 = c.rho_water*(psnowc(jk+1)/prhofirn(jk+1) + psnic(jk+1)/c.rho_ice) # layer 2
        delz = (delz1+delz2)/2             # Midpoint-to-midpoint distance
        
        dhdz = (h2-h1)/delz
        
        k1 = kF( Theta1,pdgrain(jk)  ,prhofirn(jk)  ,psnic(jk)  ,psnowc(jk), c)
        k2 = kF( Theta2,pdgrain(jk+1),prhofirn(jk+1),psnic(jk+1),psnowc(jk+1), c)
        k = (k1+k2)/2                      # Average k between the two layers

        # Initial flux according to Hirashima eqn (1)
        q0 = max(k*(dhdz + 1),0)          #
        
        # Total time-step flux according to Hirashima eqn (23)
        # and make sure it doesn't exceed available water in
        # the upper of the two layers:
        if (qlim > 0)        
            darcy_fluxes(jk) = min( pslwc(jk) , qlim*(1-exp(-q0/qlim*c.zdtime)) )
        else
            darcy_fluxes(jk) = 0
        
        
    




  # ###########################################
  # PLA Darcy (Aug 2016)
  def [liqmax] =  CLliqF(rhosin, c)
  # ###########################################

  # Calculates per pore space volume irreducible liquid water content
  # following the parameterization by Coleou and Lesaffre (1998)

    # Cap input rho_snow at pore-close off value to avoid infinite numbers
    rhos = min(rhosin,c.rho_pco)
    # Porosity
    P  = 1-rhos/c.rho_ice
    # Per mass parameterized Irreducible LWC
    Wm = 0.057 *  P/max(1-P, c.smallno) + 0.017
    # Per pore space volume
    liqmax = Wm / ( P/max(1-P, c.smallno) ) \
        * c.rho_ice/c.rho_water / max(1-Wm, c.smallno)
           
  
  def [pslwc] = hetero_percol (prhofirn, psnowc , psnic, pslwc, pdgrain, c)
# hetero_percol
#   syntax:
# [prhofirn, psnowc , psnic, pslwc, pdgrain] =\
#         hetero_percol (prhofirn, psnowc , psnic, pslwc, \
#         pdgrain, c)
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
#         pdgrain - Vector of layer grain size. see graingrowth def.
#
#         pTdeep - lower boundary temperature (K). Defined by the constant
#         T_ice_AWS so far.
#
#         zrogl - total amount of run off (mm weq)
#
#         c - Structure containing all the physical, site-depant or user
#         defined variables.
#
#   output:
#          updated [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain, zrogl]
#
#   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk)
#=========================================================================
avail_water = zeros(size(prhofirn))

# in the formulation of Marchenko et al. 2017, the liquid water is
# distributed according to a probability funcion accross the layers before
# the standard percolation scheme takes over. In theory the heterogeneous
# flow could start from any layer at depth, routing water even deeper
# That is why we commented:
# for jk = 1:c.jpgrnd-1
# and work with all the water located in the 1st layer
jk = 1

# calculating the amount of water that can be held by capillary forces
    if ( c.calc_CLliq )
         liqmaxloc = CLliqF(prhofirn(jk),c)
      else
         liqmaxloc = c.liqmax
    

    liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1)
    potret    = max( liqmaxM* psnowc(jk) , 0 )
# and what is in excess
    liqexcess = pslwc(jk) - potret
    avail_water(jk) = max ( liqexcess , 0 )
    
    #if there is available water
    if avail_water(jk)>c.smallno
        # we draw from a binomial distribution of mode hetero_percol_p to
        # know whether this layer leads to heterogeneous percolation
        if binornd(1,c.hetero_percol_p)
#             fprintf('Piping occuring\n')
            percol_water = c.hetero_percol_frac * avail_water(jk)
            #note: dflux is positive when it leaves the layer
            
            #determine random depth between current layer and maximum depth
            #range
#             thickness_weq = psnowc + snic +slwc
            thickness_act = psnowc.*(c.rho_water./prhofirn)+ psnic.*(c.rho_water/c.rho_ice)
            depth_act=zeros(size(thickness_act))
            depth_act(1)= thickness_act(1)/2 
            for i = 2:size(thickness_act,1)
                depth_act(i)=sum(thickness_act(1:i-1)) + thickness_act(i)/2 
            

            depth_current = depth_act(jk)
            depth_dest =  c.hetero_percol_dist # * rand(1,1) + depth_current
            # find index of layer at closest depth to the destination of
            # percolation
            [~, ind_dest] = min(abs(depth_act-depth_dest))
            
            sav.depth_current = depth_current
            sav.depth_dest = depth_dest
#             fprintf('from #0.2f m deep to #0.2f m\n',depth_current,depth_dest)

## In the percolation scheme of Marchenko et al. (2017) the meltwater is 
# redistributed from the surface down to a specified destination depth
# according to a probability def

# in an uniform distribution def, the subsurface layers receive water
# proportional to their thickness:
# figure
# plot(pslwc,-depth_act)
# hold on
# pslwc(jk) = pslwc(jk) - avail_water(jk)
# frac_received = thickness_act(jk:ind_dest)./sum(thickness_act(jk:ind_dest))
# pslwc(jk:ind_dest) = pslwc(jk:ind_dest) + avail_water(jk)*frac_received

# plot(pslwc,-depth_act)
# ylim([-10 0])
# pause

## an alternative is to go through the stratigraphy and stop when there is a
# gradient in grain size or an ice layer
            for ii = jk:ind_dest-1
                #here we test all the layer through which the pipe travels

                if prhofirn(ii+1)> 800
                    fprintf('ice - ')
                    break
#                 elseif psnic(ii+1)> c.ice_threshold
#                     # plus relative to frozen mass 
#                     #if there is too much ie in the next layer
#                     fprintf('ice - ')
#                     break
#                 elseif ThetaF(pslwc(ii+1), psnowc(ii+1), prhofirn(ii+1), c) >= 1
#                     #if there is a saturated layer piping can't go through
#                     #it
#                     fprintf('sat - ')
#                     break
# PERMEABILITY? h?
#                 elseif pdgrain(ii)-pdgrain(ii+1) < c.crit_diff_grain
#                     fprintf('grain - ')
#                     break
                
            
            if ind_dest ~= ii
                fprintf('Piping stopped at layer #i instead of #i\n',ii, ind_dest)
            
            ind_dest = ii
            
            jj = ind_dest
            while and(jj>jk , percol_water > c.smallno)
                # Calculate water in destination layer, when this is at saturation (Theta = 1):
                plsat = psnowc(jj)*c.rho_water/c.rho_ice *\
                    (c.rho_ice/prhofirn(jj) - 1)
                # Do not allow flux to be greater than plsat-pl in next layer.
                dflux = min (max(0, plsat-pslwc(jj)), percol_water)
                                
                pslwc(jk) = pslwc(jk) - dflux
                pslwc(jj) = pslwc(jj) + dflux
                percol_water = percol_water - dflux
                jj = jj-1
            
        
    
# 



def [hHirF] = hHirF(Theta,d, c)

#     ! The F in the name hHirF stands for "def"
# 
#     ! Calcuc.lates hydraulic suction h (in m) ac.ccording to 
#     ! Hirashima et al 2010.

alpha = 7.3*exp(1.9) #     ! Hirashima (15)
n = nHirF(d)
m=1-1/n #               ! Hirashima (9)

Theta_nozero = max(Theta,c.smallno) # ! To avoid divide-by-zero
hHirF = 1/alpha * (Theta_nozero^(-1/m)-1)^(1/n) # ! Hirashima (9)



def [kF] = kF(Theta,d,rhos,pi,ps, c)
#     ! The F in the name stands for "def"
#     ! Permeability as in Hirashima 2010, using units on Ks as in
#     ! Calonne (2012) and correction for ice layering using Colbeck 1975 (eqn 32)

#     ! Saturated hydraulic conductivity of snow (m/s) as fitted by Shimizu (1970). Here
#     ! units are as in Calonne et al (2012) and d/1000 is the grain diameter in m. It
#     ! corresponds to H10-eqn3 with units as in Calonne
# ks_shim = c.g/c.nuw*0.077*(d/1000.)^2. * exp(-7.8e-3*rhos) 
# Permeability according to Calonne (2012)
ks = 3*(d/2000)^2*c.g/c.nuw*exp(-0.013*rhos)

# fprintf('Permeability #f #f\n',ks/c.g*c.nuw/(d/2000)^2, ks_shim/c.g*c.nuw/(d/2000)^2)

# Unsaturated correction to conductivity. Hirashima eqn 10
n = nHirF(d)
m=1-1/n               # Hirashima (9)
kr = Theta^0.5 * ( 1.-(1.-Theta^(1./m))^m )^2 # Hirashima eqn 10

# Unsaturated hydraulic conductivity of the snow. Hirashima eqn 11
k11 = kr*ks

# Factor to divide k11 by to account for ice lenses within the
# layer. Colbeck 1975 eqn 32
Hs    = ps/rhos      # total depth of snow in layer
Hi    = pi/c.rho_ice    # total depth of ice in layer
fsnow = Hs / (Hs+Hi) # Fraction of layer that is snow
if (k11 > 0)
    k22factor = fsnow + (1-fsnow) * (1.+c.whwice)/(c.kice/k11 + c.whwice)
else
    # If k11 is 0 (because Theta is 0) then set to 1, so kF = 0/1 = 0 in next line
    k22factor = 1

# Effective hydraulic conductivity (in vertical direction perpicular to
# ice lenses) in snow-ice-sublayered model. Colbeck 1975 eqn 32
kF = k11/k22factor




def [nHirF] = nHirF(d)

#     ! Calcuc.lates n as in Hirashima et al 2010.
nHirF = 15.68*exp(-0.46*d)+1 # Hirashima (17)



def [qlimF] = qlimF(pl1,pl2,ps1,ps2,pi1,pi2,rhos1,rhos2,d1,d2, c)
# The F in the name stands for "def"
# Calculates Hirashima (2010) qlim in eqn 20  iteratively

# Physical layer thicknesses (m)
delz1 = c.rho_water*(ps1/rhos1 + pi1/c.rho_ice) # layer 1
delz2 = c.rho_water*(ps2/rhos2 + pi2/c.rho_ice) # layer 2
delz = (delz1+delz2)/2             # Midpoint-to-midpoint distance

# First, see if there is a solution, by moving all in upper layer down.
# If this provides positive h1-(h2+delz) then there is a solution
qlim = pl1

pl1test = pl1-qlim
pl2test = pl2+qlim
Theta1 = ThetaF(pl1test,ps1,rhos1, c)
Theta2 = ThetaF(pl2test,ps2,rhos2, c)
h1 = hHirF(Theta1,d1, c)
h2 = hHirF(Theta2,d2, c)
diff = h1 - (h2 + delz)
if ( diff < 0 )
    Conv = 1  # If moving everything isn't enough for equilibrium, move everything and do no more
    qlimOut = pl1
else
    Conv = 0  # If it is enough, we haven't converged yet and we start iteration below


if ( Conv == 0 )
    # First guess is half of water in top layer
    qlim = pl1/2
    qlimL =0
    qlimR =pl1
    
    Nmax = 11 # Number of iterations (usually 10 suffices since
    # we are halving intervals and 0.5^9 - 0.5^10 = 9.7656e-04 < 0.001)
    
    for i = 1:Nmax
        pl1test = pl1-qlim
        pl2test = pl2+qlim
        Theta1 = ThetaF(pl1test,ps1,rhos1, c)
        Theta2 = ThetaF(pl2test,ps2,rhos2, c)
        h1 = hHirF(Theta1,d1, c)
        h2 = hHirF(Theta2,d2, c)
        diff = h1 - (h2 + delz) # Difference between LHS and RHS in Hirashima eqn 20
        
        # If positive, we moved too much, and  qlim becomes new right  point
        # if negative, we moved too little, and qlim becoms new left :
        if (diff > 0 )
            qlimR = qlim
        else
            qlimL = qlim
        
        # New value is halfway between new interval  points
        qlim = (qlimR+qlimL)/2
    
    # Set final output qlim to iterated value:
    qlimOut = qlim


qlimF = qlimOut

def [ThetaF] = ThetaF(pl,ps,rhos, c)
#     ! The F in the name stands for "def"
#     ! Calculates effective water saturation.
#     ! Force Theta to be between 0 and 1
# Calculate liqmax according to Coleou and Lesaffre?
# liqmaxloc is the value of liqmax used locally here
if ( c.calc_CLliq )
  liqmaxloc = CLliqF(rhos, c)
else
  liqmaxloc = c.liqmax


if (ps > c.smallno)
    #        ! Force Theta to be between 0 and 1
    ThetaF = min(1, \
        max ( 0, \
        ( (pl / ps * rhos * c.rho_ice / c.rho_water / (c.rho_ice - rhos)\
        - liqmaxloc) / (1-liqmaxloc) ) ))
else
    #        ! If there is no snow, write out Theta=1 to avoid divide-by-zero
    ThetaF = 1



