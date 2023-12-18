# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
def error_testing():

    err(:) = 0
    o_THF = err(:)+1
    
    if c.THF_calc == 2 || c.THF_calc == 3
    	# Testing conditions required for profile method
    	# Weather variables from various origins
    	test = o_T1(:) + o_T2(:) + o_RH(:) + o_RH2(:) + o_WS(:) + o_WS2(:) ~= 0
    	err(test) = 1
    		
    	# WS below 1 m s-1
    	test = or(WS(:)<1, WS2(:)<1)
    	err(test) = 2   
    		
    	# Weather variables from same origin but T and RH measured at different height
    	test = z_T1(:) ~= z_RH(:)
    	err(test) = 3
    		
    	# WS2 <= WS
    	test = WS2(:)-WS(:)<=0
    	err(test) = 4
    		
    	if err(:) == 0       
    		[LHF(:), SHF(:), theta_2m(:), q_2m(:), ws_10m(:),Ri(:)] 			= SensLatFluxes_profile (z_T1(:),z_T2(:),			T1(:),T2(:),q(:),q2(:),WS(:),WS2(:),pres(:), c)
    
    		# Ri out of bound
    		test =  isnan(Ri(:))
    		err(test) = 5
    		
    		# Unrealistically high SHF and LHF
    		test =   abs(LHF(:))> 200 || abs(SHF(:))> 200
    		err(test) = 6
    		
    		# Surface temperature below 0K
    		test =   Tsurf(:)<=0 
    		err(test) = 7
    		
    		# Profile method output NaN (not from arleady listed errors)
    		test =   isnan(LHF(:)) ||isnan(Tsurf(:)) || isnan(SHF(:)) 
    		err(test) = 8
    		
    		LHF(err~=0)=NaN
    		SHF(err~=0)=NaN
    	
    	
    	o_THF(err==0) = 2
    	o_THF(err~=0) = 1
    	
    	if c.THF_calc == 3
    		LHF2 = LHF
    		SHF2 = SHF
    		o_THF2 = o_THF
    
    		LHF(:)=NaN
    		SHF(:)=NaN
    		o_THF(:) = 1
            

def profile_latent_heat_flux(z_r, phi_m, Ri, ro_air, c):
	# def calculating the Latent Heat Flux QE as developped by Jason Box
	# translated from IDL by Baptiste Vandecrux
	# b.vandecrux@gmail.com
	# ------------------------------------------- 
    return ro_air * c.L_vap * c.kappa**2 * z_r**2 * dq/dz * du/dz / phi_m**2


def profile_richardson(z1,z2,T1,T2,q1,q2,u1,u2,p,c):
	# def giving the richardson number developped by Jason Box
	# translated from IDL by Baptiste Vandecrux
	# b.vandecrux@gmail.com

    p0 = 1000#standard pressure for the calculation of potential temperature
    
	# The virtual temperature is the temperature that dry dry air 
	# would have if its pressure and density were equal to those 
	# of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1 = T1 * (1.+0.61*q1/1000.)# virtual temperature
    Tv2 = T2 * (1.+0.61*q2/1000.)
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

    Ri(dudz <= 0) = NaN
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
	return ro_air * c.c_pd * c.kappa**2 * z_r**2 * dtheta  / dz * du / dz  /  phi_m**2


def   profile_SensLatFluxes (z1,z2,t1,t2,q1,q2,u1,u2,p, c)
	# The virtual temperature is the temperature that dry dry air 
	# would have if its pressure and density were equal to those 
	# of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1= t1*(1.+0.61*q1)
    Tv2= t2*(1.+0.61*q1)
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
    Ri = c.g  /  Tv_avg * dtheta  /  dz * (du / dz)**(-2)

    Ri(dudz <= 0) = NaN
    if abs(Ri) > 2
        Ri = NaN
    
    
	# Stability criteria
    stab = stability(Ri)

	# air density
    ro_air = p *100 /c.R_d/ Tv_avg 

    [QH] = sensible_heat_flux (z_r,stab, Ri, ro_air,c)
    [QE]= latent_heat_flux(z1,z2,q1,q2,u1,u2, stab, Ri, ro_air, c)

    QH(Ri >= 0.2) & ( Ri < 999) = 0
    QE(Ri >= 0.2) & ( Ri < 999) = 0
    
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
    if np.sum(critical) > 0 
        Ri_temp(critical) = 0.2
    

    too_unstable = (Ri_temp <= -0.1)
    if np.sum(too_unstable) > 0
        Ri_temp(too_unstable) = -1
    

	# STABLE CASE
    phi_m = NaN(size(Ri_temp))
    stable = and(Ri_temp >= 0 , Ri_temp <= 0.2)
	# according to Box and Steffen (2001)
    phi_m(stable) = (1-5.2*Ri_temp(stable))**(-0.5)  

	# UNSTABLE CASE
    unstable = and(Ri_temp >= -1, Ri_temp < 0)
    phi_m(unstable) = (1-18*Ri_temp(unstable))**(-0.25)# Box and Steffen (2001)
    return phi_m

