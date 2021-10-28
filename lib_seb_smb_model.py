
def UpdateSnowThickness(snowthick,z_icehorizon, k, j, c):
# UpdateSnowThickness: Calculates variation of snow thickness since initial 
# conditions (snowthick_AWS).
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to python by Baptiste Vandecrux (bav@geus.dk)
#==========================================================================
	if k==1:
		#Initial snow depth
		snowthick(1,j) = c.snowthick_ini 

	#Snow to ice transition
	z_icehorizon = min(c.z_ice_max, floor(snowthick(1,j)/c.dz_ice))

	# Update BV2017: sensor height from input data
	if k>1      
		snowthick(k,j) = snowthick(k-1,j)# will be updated at  of time loop
	return snowthick, z_icehorizon  


def SurfEnergyBudget (SRnet, LRin, Tsurf, k_eff, thick_first_lay, T_ice, T_rain,    dTsurf, EB_prev, SHF, LHF, rainfall, c):
# SurfEnergyBudget: calculates the surface temperature (Tsurf) and meltflux
# from the different elements of the energy balance. The surface
# temperature is adjusted iteratively until equilibrium between fluxes is
# found.
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to python by Baptiste Vandecrux (bav@geus.dk)
#==========================================================================
    stop =0
# SURFACE ENERGY BUDGET -----------------------------------------------------------------------

    meltflux = SRnet(1) - SRnet(2) + LRin         - c.em * c.sigma * Tsurf.**4 - (1 - c.em) * LRin         + SHF + LHF         -(k_eff(1)) * (Tsurf- T_ice(2)) / thick_first_lay         + c.rho_water * c.c_w(1) * rainfall * c.dev         / c.dt_obs *( T_rain - c.T_0)
        
    if meltflux >= 0 && Tsurf == c.T_0
# stop iteration for melting surface
        stop =1     
        return
    
    
    if abs(meltflux) < c.EB_max
# stop iteration when energy components in balance
        stop =1   
        meltflux = 0
        return
    


    if meltflux/EB_prev < 0
        dTsurf = 0.5*dTsurf 
# make surface temperature step smaller when it overshoots EB=0
    
    EB_prev = meltflux

#Update BV
    if meltflux < 0
        Tsurf = Tsurf - dTsurf 
    else          
        Tsurf = min(c.T_0,Tsurf + dTsurf)
    return meltflux, Tsurf, dTsurf, EB_prev, stop



def SRbalance (SRout_mdl, SRin, SRnet, z_icehorizon,snowthick, T_ice, rho, j, k, c):

# SRbalance: Calculates the amount of Shortwave Radiation that is 
# penetrating at each layer (SRnet).
# uses it to warm each layer and eventually calculates the melt that is
# produced by this warming
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to python by Baptiste Vandecrux (bav@geus.dk)
#==========================================================================
#extinction coefficient of ice 0.6 to 1.5 m-1
#extinction coefficient of snow 4 to 40 m-1
# Greufell and Maykut, 1977, Journal of Glacionp.logy, Vol. 18, No. 80, 1977

#radiation absorption in snow
# SRnet(snow_layer) = (SRin - SRout_mdl)#     *exp(-ext_snow*depth(snow_layer))
# 
#radiation absorption in ice layers underneath the snowpack
#  SRnet(ice_layer) = (SRin-SRout).*#         exp(-ext_snow*snowthick).*#         exp(-ext_ice*(depth(ice_layer) - snowthick))

# NOT USED AT THE MOMENT
# should be updated for non regular layer thickness
    
if c.elev_bins ~= 1 && k > 1
# in a transect, SRout is calculated using SRin (calculated from SRin_AWS)
# and fixed values for albedo for snow and ice
    if snowthick(k,j) > 0
        SRout_mdl(k,j) = c.alb_snow*SRin(k,j) 
    else
        SRout_mdl(k,j) = c.alb_ice*SRin(k,j)
    

#radiation absorption in snow

SRnet(1:(z_icehorizon+1)) = (SRin(k,j) - SRout_mdl(k,j))    *exp(-c.ext_snow*(0:z_icehorizon)*c.dz_ice)

#radiation absorption in underlying ice
if z_icehorizon < c.z_ice_max
    SRnet((z_icehorizon+2):c.z_ice_max+1) = (SRin(k,j)-SRout_mdl(k,j)).*        exp(-c.ext_snow*snowthick(k,j)).*        exp(-c.ext_ice*(((z_icehorizon+2):(c.z_ice_max+1))*        c.dz_ice - snowthick(k,j)))

# specific heat of ice 
#(perhaps a slight overestimation for near-melt ice temperatures (max 48 J/kg/K))
if k==1
    c_i = 152.456 + 7.122 * T_ice(:,1,j)    
#snow & ice temperature rise due to shortwave radiation absorption
else
# Specific heat of ice (a slight overestimation for near-melt T (max 48 J kg-1 K-1))
    c_i = 152.456 + 7.122 * T_ice(:,k-1,j)
    
    T_ice(1:c.z_ice_max,k,j) = T_ice(1:c.z_ice_max,k-1,j) +        c.dt_obs./c.dev./rho(1:c.z_ice_max,k)./c_i(1:c.z_ice_max).*        (SRnet(1:c.z_ice_max)-SRnet(2:(c.z_ice_max+1)))./c.dz_ice
    
#      a=#         (SRnet(1:c.z_ice_max)-SRnet(2:(c.z_ice_max+1)))./c.dz_ice
    
    T_ice(c.z_ice_max+1,k,j) = T_ice(c.z_ice_max+1,1,j)



# finding where/how much melt occurs
subsurfmelt = (T_ice(:,k,j) > c.T_0)
nosubsurfmelt = (T_ice(:,k,j) <= c.T_0)
dT_ice = T_ice(:,k,j) - c.T_0
dT_ice(nosubsurfmelt) = 0

meltflux_internal_temp = rho(:,k) .* c_i .* dT_ice / c.dt_obs * c.dev * c.dz_ice
meltflux_internal = sum(meltflux_internal_temp(1:c.z_ice_max))

dH_melt_internal = zeros (c.jpgrnd,1)
dH_melt_internal(:) = -c_i.*dT_ice./c.L_fus * c.dz_ice
dH_melt_internal(1) = 0

if sum(subsurfmelt) > 0
    T_ice(subsurfmelt,k,j) = c.T_0# removing non-freezing temperatures

# Reduce sub-surface density due to melt? Will most likely cause model instability\      

return SRout_mdl, SRnet, T_ice, meltflux_internal, dH_melt_internal


def RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c):

u_star = c.kappa * WS / (np.log(z_WS/z_0) - psi_m2 + psi_m1)  

# rough surfaces: Smeets & Van den Broeke 2008
Re = u_star * z_0 / nu       

z_h = z_0 * exp(1.5 - 0.2*np.log(Re) - 0.11*(np.log(Re))**2) 

if z_h < 1e-6
     z_h = 1e-6

z_q = z_h

return z_h, z_q, u_star, Re

def  SmoothSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c):

u_star = c.kappa*WS/(np.log(z_WS/z_0)- psi_m2 + psi_m1)  
Re = u_star * z_0 / nu 

if Re <= 0.135
    range = 1
elseif and(Re > 0.135 , Re < 2.5)
    range = 2
elseif Re >= 2.5
    range = 3
else
    disp('ERROR')
    disp(Re)
    disp(u_star)
    disp(z_WS)
    disp(psi_m2)
    disp(psi_m1)
    disp(nu)
    disp(range)


# smooth surfaces: Andreas 1987
z_h =  z_0 * exp(c.ch1(range) + c.ch2(range)*np.log(Re) + c.ch3(range)*(np.log(Re))**2)  
z_q = z_0 * exp(c.cq1(range) + c.cq2(range)*np.log(Re) + c.cq3(range)*(np.log(Re))**2)

    if z_h < 1e-6
         z_h = 1e-6
    
    if z_q < 1e-6
         z_q = 1e-6
    
return z_h, z_q, u_star, Re


def SensLatFluxes_bulk (WS, nu, q, snowthick, Tsurf,theta, theta_v , pres,rho_atm,  z_WS, z_T, z_RH, z_0, c):
# SensLatFluxes: Calculates the Sensible Heat Flux (SHF), Latent Heat Fluxes
# (LHF) and Monin-Obhukov length (L). Corrects for atmospheric stability.
#
# see Van As et al., 2005, The Summer Surface Energy Balance of the High 
# Antarctic Plateau, Boundary-Layer Meteoronp.logy.
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to python by Baptiste Vandecrux (bav@geus.dk)
#==========================================================================
psi_m1 = 0
psi_m2 = 0

# will be updated later
theta_2m = theta
q_2m     = q
ws_10m    =  WS
        
if WS > c.WS_lim
# Roughness length scales for snow or ice - initial guess
    if WS<c.smallno
        z_h = 1e-10
        z_q = 1e-10
    else
        if snowthick > 0
            [z_h, z_q, u_star, Re] = SmoothSurf(WS,z_0, psi_m1, psi_m2, nu, z_WS, c)
        else
            [z_h, z_q, u_star, Re] = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
        
    


    es_ice_surf = 10.**(-9.09718 * (c.T_0 / Tsurf - 1.)         - 3.56654 * np.log10(c.T_0 / Tsurf) + 0.876793 * (1. - Tsurf / c.T_0)         + np.log10(c.es_0))
    q_surf  = c.es * es_ice_surf/(pres-(1-c.es)*es_ice_surf)
    L = 10e4

    if theta >= Tsurf && WS >= c.WS_lim# stable stratification
# correction from Holtslag, A. A. M. and De Bruin, H. A. R.: 1988, ‘Applied Modelling of the Night-Time
# Surface Energy Balance over Land’, J. Appl. Meteorol. 27, 689–704.
        for i=1:c.iter_max_flux
            psi_m1 = -(c.aa*z_0/L  +  c.bb*(z_0/L-c.cc/c.dd)*exp(-c.dd*z_0/L)  + c.bb*c.cc/c.dd)
            psi_m2 = -(c.aa*z_WS/L + c.bb*(z_WS/L-c.cc/c.dd)*exp(-c.dd*z_WS/L) + c.bb*c.cc/c.dd)
            psi_h1 = -(c.aa*z_h/L  +  c.bb*(z_h/L-c.cc/c.dd)*exp(-c.dd*z_h/L)  + c.bb*c.cc/c.dd)
            psi_h2 = -(c.aa*z_T/L  +  c.bb*(z_T/L-c.cc/c.dd)*exp(-c.dd*z_T/L)  + c.bb*c.cc/c.dd)
            psi_q1 = -(c.aa*z_q/L  +  c.bb*(z_q/L-c.cc/c.dd)*exp(-c.dd*z_q/L)  + c.bb*c.cc/c.dd)
            psi_q2 = -(c.aa*z_RH/L  +  c.bb*(z_RH/L-c.cc/c.dd)*exp(-c.dd*z_RH/L)  + c.bb*c.cc/c.dd)
#     if psi_m2 <-10000
#         df = 0
#     
            if WS<c.smallno
                z_h = 1e-10
                z_q = 1e-10
            else              
                if snowthick > 0
                    [z_h, z_q, u_star, Re] = SmoothSurf(WS,z_0, psi_m1, psi_m2, nu, z_WS, c)
                else
                    [z_h, z_q, u_star, Re] = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
                
            

            th_star =  c.kappa * (theta - Tsurf) / (np.log(z_T / z_h) - psi_h2 + psi_h1) 
            q_star  =  c.kappa * (q - q_surf) / (np.log(z_RH / z_q) - psi_q2 + psi_q1) 
            SHF  = rho_atm * c.c_pd  * u_star * th_star
            LHF  = rho_atm * c.L_sub * u_star * q_star

            L_prev  = L
#             L    = u_star**2 * theta_v  / ( 3.9280 * th_star*(1 + 0.6077*q_star))
        L    = u_star**2 * theta*(1 + ((1 - c.es)/c.es)*q)  / (c.g * c.kappa * th_star*(1 + ((1 - c.es)/c.es)*q_star))

            if L == 0 || (abs((L_prev - L)) < c.L_dif)
# calculating 2m temperature, humidity and wind speed
                theta_2m = Tsurf + th_star/c.kappa * (np.log(2/z_h) - psi_h2 + psi_h1)
                q_2m     = q_surf + q_star /c.kappa * (np.log(2/z_q) - psi_q2 + psi_q1)
                ws_10m    =          u_star /c.kappa * (np.log(10/z_0) - psi_m2 + psi_m1)
                break
            
        

    

    if theta < Tsurf && WS >= c.WS_lim# unstable stratification
# correction defs as in 
# Dyer, A. J.: 1974, ‘A Review of Flux-Profile Relationships’, Boundary-Layer Meteorol. 7, 363– 372.
# Paulson, C. A.: 1970, ‘The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer’, J. Appl. Meteorol. 9, 857–861.

        for i=1:c.iter_max_flux
            x1      = (1 - c.gamma * z_0  / L)**0.25
            x2      = (1 - c.gamma * z_WS / L)**0.25
            y1      = (1 - c.gamma * z_h  / L)**0.5
            y2      = (1 - c.gamma * z_T  / L)**0.5
            yq1     = (1 - c.gamma * z_q  / L)**0.5
            yq2     = (1 - c.gamma * z_RH  / L)**0.5
            psi_m1  = np.log( ((1 + x1)/2)**2 * (1 + x1**2)/2) - 2*atan(x1) + pi/2        
            psi_m2  = np.log( ((1 + x2)/2)**2 * (1 + x2**2)/2) - 2*atan(x2) + pi/2
            psi_h1  = np.log( ((1 + y1)/2)**2 )
            psi_h2  = np.log( ((1 + y2)/2)**2 )
            psi_q1  = np.log( ((1 + yq1)/2)**2 ) 
            psi_q2  = np.log( ((1 + yq2)/2)**2 )

            if WS<c.smallno
                z_h = 1e-10
                z_q = 1e-10
            else
                if snowthick > 0
                    [z_h, z_q, u_star, Re] = SmoothSurf(WS,z_0, psi_m1, psi_m2, nu, z_WS, c)
                else
                    [z_h, z_q, u_star, Re] = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
                
            

            th_star = single (c.kappa * (theta - Tsurf) / (np.log(z_T / z_h) - psi_h2 + psi_h1))
            q_star  = single ( c.kappa * (q  - q_surf) / (np.log(z_RH / z_q) - psi_q2 + psi_q1) )
            SHF  = single ( rho_atm * c.c_pd  * u_star * th_star)
            LHF  = single( rho_atm * c.L_sub * u_star * q_star)

            L_prev  = L
            L    = u_star**2 * theta*(1 + ((1 - c.es)/c.es)*q) /                (c.g * c.kappa * th_star*(1 + ((1 - c.es)/c.es)*q_star))
# if i ==1
#     figure
#     title('Unstable')
# 
#     hold on
# 
# scatter(i,L)
            if abs((L_prev - L)) < c.L_dif
# calculating 2m temperature, humidity and wind speed
                theta_2m = Tsurf + th_star/c.kappa * (np.log(2/z_h) - psi_h2 + psi_h1)
                q_2m     = q_surf + q_star /c.kappa * (np.log(2/z_q) - psi_q2 + psi_q1)
                ws_10m    =          u_star /c.kappa * (np.log(10/z_0) - psi_m2 + psi_m1)
               
                break
            

        
    

else
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
    psi_q1 = 0
    psi_q2 = -999
    
# calculating 2m temperature, humidity and wind speed
    theta_2m = theta
    q_2m     = q
    ws_10m    =  WS 

def L, LHF, SHF, theta_2m, q_2m , ws_10m, Re


def SpecHumSat(RH, T, pres, c):

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
#==========================================================================


# SPECIFIC HUMIDITY & SATURATION -----------------------------------------------------------------------
es_wtr = 10.**(-7.90298*(c.T_100./T-1)     + 5.02808 * np.log10(c.T_100./T) \# saturation vapour pressure above 0 C (hPa)
    - 1.3816E-7 * (10.**(11.344*(1.-T./c.T_100))-1.)     + 8.1328E-3*(10.**(-3.49149*(c.T_100./T-1)) -1.) + np.log10(c.es_100))

es_ice = 10.**(-9.09718 * (c.T_0 ./ T - 1.)     - 3.56654 * np.log10(c.T_0 ./ T) +     0.876793 * (1. - T ./ c.T_0) + np.log10(c.es_0))# saturation vapour pressure below 0 C (hPa)
q_sat = c.es * es_wtr./(pres-(1-c.es)*es_wtr)# specific humidity at saturation (incorrect below melting point)

freezing = find(T < c.T_0)# replacing saturation specific humidity values below melting point
if sum(freezing) > 0
    q_sat(freezing) = c.es * es_ice(freezing)./(pres(freezing)-(1-c.es)*es_ice(freezing))

# supersaturated = find(RH > 100)# replacing values of supersaturation by saturation
# if sum(supersaturated) > 0
#     RH(supersaturated) = 100
# 
q = RH.*q_sat/100# specific humidity in kg/kg
return RH, q

def profile_latent_heat_flux(z_r, phi_m, Ri, ro_air, c):
# def calculating the Latent Heat Flux QE as developped by Jason Box
# translated from IDL by Baptiste Vandecrux
# b.vandecrux@gmail.com
# ------------------------------------------- 
    return ro_air * c.L_vap * c.kappa**2 * z_r**2 * dq/dz * du/dz / phi_m.**2


def profile_richardson(z1,z2,T1,T2,q1,q2,u1,u2,p,c):
# def giving the richardson number developped by Jason Box
# translated from IDL by Baptiste Vandecrux
# b.vandecrux@gmail.com

    p0 = 1000#standard pressure for the calculation of potential temperature
    
# The virtual temperature is the temperature that dry dry air 
# would have if its pressure and density were equal to those 
# of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1 = T1 .* (1.+0.61*q1/1000.)# virtual temperature
    Tv2 = T2 .* (1.+0.61*q2/1000.)
    Tv_avg = (Tv1+Tv2)/2

# The potential temperature  is the temperature that an unsaturated
# parcel of dry air would have if brought adiabatically and 
# reversibly from its initial state to a standard pressure, p0,
# typically 1000 hPa
    VPT1 = Tv1*(p0/p)**c.kappa_poisson#virtual potential temperatures
    VPT2 = Tv2*(p0/p)**c.kappa_poisson

    DVPT=VPT2-VPT1
    du = u2-u1
    dz = z2-z1
    
    Ri = c.g / Tv_avg * DVPT / dz * (du/dz)**(-2)

    Ri(dudz <= 0 ) = NaN
    if abs(Ri) > 2
        Ri = NaN
	return Ri
    


def profile_sensible_heat_flux (z_r,phi_m, Ri, ro_air,c):
# def calculating the Sensible Heat Flux QH usingthe profile method
# as developped by Jason Box.
# translated from IDL by Baptiste Vandecrux
# b.vandecrux@gmail.com
# ------------------------------------------- 

# phi_h needs to be tuned up for very unstable conditions according to Box and Steffen (2001)
phi_m(Ri<-0.03) = phi_m(Ri<-0.03)*1.3
return ro_air .* c.c_pd .* c.kappa**2 .* z_r.**2 .* dtheta ./dz .* du./dz ./ phi_m.**2


def   profile_SensLatFluxes (z1,z2,t1,t2,q1,q2,u1,u2,p, c)
        
# The virtual temperature is the temperature that dry dry air 
# would have if its pressure and density were equal to those 
# of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1=  t1*(1.+0.61*q1)
    Tv2=  t2*(1.+0.61*q1)
    Tv_avg=(Tv1 + Tv2)/2.

# The potential temperature  is the temperature that an unsaturated
# parcel of dry air would have if brought adiabatically and 
# reversibly from its initial state to a standard pressure, p0,
# typically 1000 hPa
# Here we calculate the virtual potential temperature (dry)
    p0 = 1000#standard pressure for the calculation of potential temperature

    theta_v1=t_v1*(p0/p)**c.kappa_poisson
    theta_v2=t_v2*(p0/p)**c.kappa_poisson
    dtheta = theta_v2-theta_v1
    du = u2-u1
    dz = z2-z1
    z_r = dz/np.log(z2/z1)
# z_r = sqrt(z1*z2)

# Richardson number
    Ri = c.g ./ Tv_avg .* dtheta ./ dz .* (du./dz)**(-2)

    Ri(dudz <= 0 ) = NaN
    if abs(Ri) > 2
        Ri = NaN
    
    
# Stability criteria
    stab = stability(Ri)

# air density
    ro_air =  p *100 /c.R_d/ Tv_avg 

    [QH] = sensible_heat_flux (z_r,stab, Ri, ro_air,c)
    [QE]=  latent_heat_flux(z1,z2,q1,q2,u1,u2, stab, Ri, ro_air, c)

    QH(Ri >= 0.2 && Ri < 999) = 0
    QE(Ri >= 0.2 && Ri < 999) = 0
    
    T2m = t2
    q_2m = q2
    ws_10m = u2
	return QE,QH, T2m, q_2m, ws_10m, Ri
    
	
def profile_stability(Ri):
# def calculating the correction parameters for the similarity 
# developped by Jason Box
# translated from IDL by Baptiste Vandecrux
# b.vandecrux@gmail.com

    Ri_temp=Ri

    critical = and(Ri_temp >= 0.2 , Ri_temp < 999)
    if sum(critical) > 0 
        Ri_temp(critical) = 0.2
    

    too_unstable = (Ri_temp <= -0.1)
    if sum(too_unstable) > 0
        Ri_temp(too_unstable) = -1
    

# -------------------------------- STABLE CASE
    phi_m = NaN(size(Ri_temp))
    stable = and(Ri_temp >= 0 , Ri_temp <= 0.2)
# according to Box and Steffen (2001)
    phi_m(stable) = (1-5.2*Ri_temp(stable))**(-0.5)  

# -------------------------------- UNSTABLE CASE
    unstable = and(Ri_temp >= -1, Ri_temp < 0)
    phi_m(unstable) = (1-18*Ri_temp(unstable))**(-0.25)# Box and Steffen (2001)
	return phi_m


