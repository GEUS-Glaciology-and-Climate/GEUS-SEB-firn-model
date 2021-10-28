def [rho_snow] = IniRhoSnow(T, WS, elev, c)
#IniRhoSnow: Initialization of fresh snow density using different
#parametreization. If modified, check that the chosen parametrization is
#sent properly to the subsurface scheme.
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
#==========================================================================

# INITIALIZATION OF SNOW DENSITY 
rho_snow = zeros(size(T))
mean_T = zeros(c.elev_bins)
for j=1:c.elev_bins
    mean_T = mean(T(:,j))


for j=1:c.elev_bins
#     mean_T(j) = mean(T(j,:))
    switch c.rho_snow_scheme
        case 0
            rho_snow(:,j) = c.fresh_snow_dens #constant in time and place
        case 1
            
            rho_snow(:,j) = 625. + 18.7*(mean_T(j)-c.T_0) + 0.293*(mean_T(j)-c.T_0)^2
            # surface snow density by Reeh et al. 2005
        case 2
            rho_snow(:,j) = 481. + 4.834*(mean_T(j)-c.T_0)
            # surface snow density by Kuipers-Munneke et al. 2015
        case 3
            rho_snow(:,j) = 369.58 -0.02298*elev(j) 
            # PLA RFO fresh snow density (May 2015): BV2017 using
            # elevation-based parametrization
        case 4
            rho_snow(:,j) = 350. #default value everywhere
            ind = (T(:,j) >= 258.16)
            rho_snow(ind,j) = 50 + 1.7 * (T(ind,j) - 258.16).^1.5
            clear ind
            ind = and((T(:,j) >= 258.16) , (WS(:,j) >= 5))
            rho_snow(ind,j) = 50 + 1.7 * (T(ind,j) - 258.16).^1.5 +\
                25 + 250 *(1-exp(-0.2*(WS(ind) - 5)))
            # surface snow density by Liston et al., 2007
            # only works above -15 degC
    



# This script helps to generate a density profile to be used as initial
# state for the subsurface model
#
# Author: Baptiste Vandecrux (bava@byg.dtu.dk)
# ========================================================================

# clear all
# close all
load ../matlab_defs/Core_all


station = 'KAN_U'
# Plot in the site you want
# PlotCore(Core,'CoreNum',1)

# Find the index of the core you want to use
# CoreList(CoreAvg)
# search by name rather than by index since index can change in the future

switch station
    case 'CP1'
        ind = FindCore(Core,'Name','CORE 6945')
    case 'DYE-2'
        ind = FindCore(Core,'Name','DYE2 1998 core B')
    case 'DYE-2_HQ'
        ind = FindCore(Core,'Name','core_10_2016')
    case 'KAN-U'
        ind = FindCore(Core,'Name','core_1_2012')
    case 'NASA-SE'
        ind = FindCore(Core,'Name','CORE 6642 (B)')
    case 'EGP'
        ind = FindCore(Core,'Name','NEGIS')
#        ind_new = length(Core)+1
#        Core{ind_new} = Core{ind1}
#        Core{ind_new}.Info.Name = 'NASA-SE_bapt'
#        Core{ind_new}.Data.Density(50:110) = NaN
#        ind_nan = isnan(Core{ind_new}.Data.Density)
#        Core{ind_new}.Data.Density(ind_nan) = \
#            interp1(Core{ind_new}.Data.Depth(~ind_nan),\
#            Core{ind_new}.Data.Density(~ind_nan),\
#            Core{ind_new}.Data.Depth(ind_nan))
#        figure
#        OverlapPlot(Core,[ind1 ind_new])
#         ind = ind_new
        
#         ind2=FindCore(Core,'NearestCodeLocation','NASA-SE')
#         figure
#         PlotCore(Core,'CoreNumber',ind2)
#        figure
#        OverlapPlot(Core,[114 23])
    case 'Summit'
        ind1 = FindCore(Core,'Name','T99_1990')
       ind2 = FindCore(Core,'Name','Grip1991Shallow')
       ind_new = length(Core)+1
       Core{ind_new} = Core{ind2}
       Core{ind_new}.Info.Name = 'T99+GRIP_1990'
       Core{ind_new}.Data.Density(1:length(Core{ind1}.Data.Density)) = \
           Core{ind1}.Data.Density
       ind_nan = isnan(Core{ind_new}.Data.Density)
       Core{ind_new}.Data.Density(ind_nan) = \
           interp1(Core{ind_new}.Data.Depth(~ind_nan),\
           Core{ind_new}.Data.Density(~ind_nan),Core{ind_new}.Data.Depth(ind_nan))
       
        ind = ind_new
    case 'Miege'
        ind = 120 #or 8
    case 'NASA-U'
        ind = FindCore(Core,'Name','CORE 7347')
#         ind2 = FindCore(Core,'Name','NASA-U_Herson')

    case 'SouthDome'
        ind = FindCore(Core,'Name','S. Dome Core A')
#         ind = FindCore(Core,'Name','S. Dome Core B')
    case 'NASA-E'
        ind1 = FindCore(Core,'Name','NASA East Core A')
        ind2 = FindCore(Core,'Name','NASA East Core B')
        ind_new = length(Core) +1
        Core{ind_new} = Core{ind1}
        Core{ind_new}.Info.Name = 'NASA East combined'
        tmp = length(Core{ind2}.Data.Density)
        Core{ind_new}.Data.Density(1:tmp) = Core{ind2}.Data.Density

    case 'Saddle'
        ind = FindCore(Core,'Name','N. Dye 3 (Saddle) - B')         
        ind2 = FindCore(Core,'Name','N. Dye 3 (Saddle) - A')
        figure
        OverlapPlot(Core,[ind ind2])
    case 'TUNU-N'
        ##
#         ind2 = FindCore(Core,'Name','B19_NGT19_1994')
#         ind2 = FindCore(Core,'Name','Tunu-S7.5')
#         tunus7.5 lower than tunu1
#         B19 lower than tunu1
#         tunuN25 lower than tunu1
#         tunuE25 lower than tunu1
#         tunuW25 lower than tunu1


        ind1 = FindCore(Core,'Name','Tunu-W25')
        ind2 = FindCore(Core,'Name','Tunu-E50')

        ind = FindCore(Core,'Name','Tunu-1')

        figure
        OverlapPlot(Core,[ind ind2])
        ##
    case 'NGRIP'
        ind = FindCore(Core,'Name','NG97S2~1-3bag')
#         ind = FindCore(Core,'Name','NGRIP2001S5')
    case 'GITS'
        ind = FindCore(Core,'Name','Camp Century')
#         ind = FindCore(Core,'Name','NGRIP2001S5')
    case 'KAN_U'
        ind = FindCore(Core,'Name','core_1_2013')
#         ind = FindCore(Core,'Name','NGRIP2001S5')
       


# Creating density profile and writing it into the Input folder
depth = Core{ind}.Data.Depth/100
density = Core{ind}.Data.Density
ice_perc = Core{ind}.Data.Type_perc

# if strcmp(station,'KAN-U')
#     depth = depth'
#     density = density'
#     ice_perc = ice_perc'
#     ice_perc(length(ice_perc):length(density)) = 0
# 
# DensProfile = [depth, density, ice_perc]
filename = sprintf('./Input/Initial state/DensityProfile_#s_#i.csv',station,Core{ind}.Info.DateCored.Year)
    M  = [depth, density]
    M_table = array2table(M,'VariableName', {'depth_m', 'density_kgm3'})

    writetable(M_table,filename,'Delimiter','')

fprintf('Initial density profile was generated from core #s and placed in Input folder.\n',Core{8}.Info.Name)
PlotCore(Core,'CoreNumber',ind)






def [T_ice, rhofirn_ini,rho,snic_ini, snowc_ini, \
    slwc_ini, graind_out] = \
    InitializationSubsurface(T_obs, depth_thermistor, T_ice, \
    time, Tsurf_ini, j, c)

# InitializationSubsurface: Sets the initial state of the sub surface parameters:
# - snow depth
# - temperature profile
# - density profile
#
#
# Author: Baptiste Vandecrux (bava@byg.dtu.dk)
#==========================================================================


## ========== Initial density profile ==================================
# Here we read the density profile from a core defined in a csv file in
# Input folder, average it on the weq depth scale used by the model, and
# give it as initial density profile for the model.  Below the depth of the
# given core, it is assumed to find ice.
#
# See script "InitialDensityProfileGenerator.m" in lib for creation of the
# csv file.
if c.retmip == 0
    switch c.station
        case {'DYE-2','DYE-2_long'}
            filename = './Input/Initial state/density/RetMIP_density_Dye-2 1998.csv'
        case {'CP1', 'Summit', 'NASA-SE', 'NASA-E', 'NASA-U', 'SouthDome', 'Saddle', 'TUNU-N'}
            filename = ['./Input/Initial state/density/DensityProfile_', c.station,'_1998.csv']
#             filename = './Input/Initial state/density/RetMIP_density_Summit 1990.csv'
        case 'NGRIP'
            filename = './Input/Initial state/density/DensityProfile_NGRIP_1997.csv'  
        case 'GITS'
            filename = './Input/Initial state/density/DensityProfile_GITS_1963.csv'  
        case 'DYE-2_HQ'
            filename = './Input/Initial state/density/RetMIP_density_Dye-2 2016.csv'
        case 'NUK_K'
            filename = './Input/Initial state/density/DensityProfile_NUK_K.csv'
        case 'EGP'
            filename = './Input/Initial state/density/DensityProfile_EGP_2012.csv'
        case 'KAN-U'
    #         if c.year(1) == 2012
                filename = './Input/Initial state/density/RetMIP_density_KAN-U 2012.csv'
    #         else
    #             filename = './Input/Initial state/density/DensityProfile_KAN-U_1989.csv'
    #         
       case 'Miege'
            filename = './Input/Initial state/density/RetMIP_density_Firn Aquifer.csv'
        otherwise
            if c.elev_AWS < 1400
                filename = './Input/Initial state/density/DensityProfile_NUK_K.csv'
            elseif c.elev_AWS < 1800
                filename = './Input/Initial state/density/RetMIP_density_Dye-2 2016.csv'
            else
                filename = './Input/Initial state/density/DensityProfile_Summit_1998.csv'
            
            warning('Missing initial density profile for requested station.')
    
else
    filename = ['.\RetMIP\Input files\density\RetMIP_density_', c.station,'.csv']

disp('File used for firn density initialization:')
disp(filename)

delimiter = ''
startRow = 1
formatSpec = '#f#f#[^\n\r]'

fileID = fopen(filename,'r')
try dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,\
    'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false)
catch me
    startRow = 2
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,\
        'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false)

fclose(fileID)

out = [dataArray{1:-1}]

depth = out(:,1) #real depth scale
olddensity = out(:,2) #density in kg/m^3 on the old real depth scale
if size(out,2)>2
    oldstrat = out(:,3) #density in kg/m^3 on the old real depth scale
else
    oldstrat=zeros(size(olddensity))


if c.perturb_rho_ini ~= 0
    old_thickness = depth
    old_thickness(2: ) = depth(2:) - depth(1:-1)
    [~, ind_20] = min(abs(depth-20))

    olddensity_avg = sum(olddensity(1:ind_20).*old_thickness(1:ind_20)) / depth(ind_20)
    perturb_olddensity = olddensity(1:ind_20)\
        + c.perturb_rho_ini * depth(ind_20)/length(1:ind_20)./ old_thickness(1:ind_20)
    perturb_olddensity_avg = sum(perturb_olddensity.*old_thickness(1:ind_20))\
        / depth(ind_20)
    perturb_olddensity(perturb_olddensity>917)=917

    f = figure('Visible','on','OuterPosition',[0 0 8 18])
    plot(olddensity,-depth)
    hold on
    plot(perturb_olddensity,-depth(1:ind_20))
    ylabel('Depth (m)')
        xlabel('Density (kg m^{-3})','Interpreter','tex')
    leg('Original','After perturbation')
    title(sprintf('Initial density perturbated\non average by #0.1f kg/m-3 \nin the upper 20 m.',perturb_olddensity_avg - olddensity_avg))

    olddensity(1:ind_20) = perturb_olddensity

#the depthweq of each layer dep on the thickness_weq of all overlying
#layers. It is therefore necessary to fill missing data. Missing density
#are thus replaced by 350 kg/m^3.
olddensity(isnan(olddensity))=350

thickweq = zeros(size(depth)) #thickness of each old layer in water eq
#calculating weq thickness of each layer
# tick_old_act =  depth
# tick_old_act(2:)  = depth(2:) - depth(1:-1)
# thickweq = thick_old_act .* olddensity /c.rho_water
thickweq(1:length(depth)-1) = (depth(2:length(depth))-depth(1:length(depth)-1)) \
    .*olddensity(1:length(depth)-1)/c.rho_water
thickweq() = (depth()-depth(-1))*olddensity()/c.rho_water

oldscale_weq = zeros(size(depth)) #depthweq of each layer in original scale
# calculating the old weq depth scale
oldscale_weq(:) = cumsum(thickweq(:))

# calculates the new depth scale in mweq
depth_weq = zeros(size(c.cdel))
depth_weq = cumsum(c.cdel)
rhofirn_ini = zeros(size(c.cdel))
snowc_ini = zeros(size(c.cdel))
snic_ini = zeros(size(c.cdel))
slwc_ini = zeros(size(c.cdel))

# we limit the new scale to the values that are given within the oldscale
newscale_weq = depth_weq(depth_weq<oldscale_weq())
olddensity = olddensity(oldscale_weq<=newscale_weq())
depth = depth(oldscale_weq<=newscale_weq())
oldstrat = oldstrat(oldscale_weq<=newscale_weq())
oldscale_weq = oldscale_weq(oldscale_weq<=newscale_weq())
newscale_weq = depth_weq(depth_weq<oldscale_weq())

# Since the oldscale is denser than the new one, we need to make the union
# of both scales and then average the density for each section of the new
# scale.
mergedscale = union(oldscale_weq,newscale_weq)
density_mergedscale = interp1(oldscale_weq,olddensity,mergedscale,'next','extrap')
ice_perc_mergedscale = interp1(oldscale_weq,oldstrat,mergedscale,'next','extrap')

# In the coming section we further process the core to increase resolution.
# The stratigraphy is used to know the percentage of ice within each core
# section. Then it is used to initiate snic. The firn remaining in each
# section is used to initiate snowc and the weight of the section (once
# subtracted the weight of ice content) is used to initiate prhofirn.

    # calculating the thickness of each layer in merged scale
    thick_mergedscale_weq = mergedscale
    thick_mergedscale_weq(2:length(mergedscale)) = mergedscale(2:) - mergedscale(1:-1)
    thick_mergedscale_m = thick_mergedscale_weq * c.rho_water \
        ./ density_mergedscale

    # calculating the volume and mass of ice contained in each section
    icevol_mergedscale_m = thick_mergedscale_m .* ice_perc_mergedscale/100
    icemass_mergedscale_kg = icevol_mergedscale_m * 800 
    # in line above 600 is taken for density of ice, it is low but it was
    # noticed that putting higher values  leads to very low density for the
    # firn left in the section.

    # calculating the volume and mass of firn left in the section excluding
    # the ice content
    firnvol_mergedscale_m = thick_mergedscale_m .* (1-ice_perc_mergedscale/100)
    firnmass_mergedscale_kg = density_mergedscale .* thick_mergedscale_m - icemass_mergedscale_kg
    firndensity_mergedscale = firnmass_mergedscale_kg ./ firnvol_mergedscale_m
    firndensity_mergedscale(firnvol_mergedscale_m==0)=500 # giving an aritificial density of 500 to the firn in sections only made of ice
    
    ind_bins = discretize(mergedscale,[0 depth_weq],'IncludedEdge','right')
    nbins = length(unique(ind_bins))
    mean_firndensity = zeros(nbins, 1)

    # now we either sum or average the ice content, snow content and firn
    # density on the model depth scale
    for ii = 1:nbins
        ind_in_bin = (ind_bins == ii)

        mean_firndensity(ii) = sum(density_mergedscale(ind_in_bin).*thick_mergedscale_m(ind_in_bin)) \
            /sum(thick_mergedscale_m(ind_in_bin))

        snowc_ini = c.cdel
        snic_ini(ii) = 0
        slwc_ini(ii) = 0
        if snowc_ini(ii)<0
            eweojjnjv
        
    
 
# giving the observed density (on appropriate scale) as initial value
# for the subsurface density profile
rhofirn_ini(1:length(mean_firndensity)) = min(max(300, mean_firndensity),900)

# fills the rest of the density profile with ice
if strcmp(c.station,'NUK_K')
        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = \
            860
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1)
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0
else
    if length(mean_firndensity)<length(c.cdel)
        x =[depth_weq(1:length(mean_firndensity)) 30 70]
        y = [mean_firndensity 830830]
        p = polyfit(x,y,2)
        fo = @(x) p(1)*x.^2 + p(2)*x + p(3)
        
        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = \
            fo(depth_weq(length(mean_firndensity)+1:length(c.cdel)))
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1) 
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0
    


rho =  (c.cdel *c.rho_water)./\
    ((snowc_ini *c.rho_water)./rhofirn_ini + (snic_ini *c.rho_water)./c.rho_ice)

# Removing densities greater than ice densities (just in case)
toodense = (rhofirn_ini > c.rho_ice)
rhofirn_ini(toodense) = c.rho_ice

# calculates the new depth scale in real m
depth_act = cumsum(c.cdel .*c.rho_water ./rhofirn_ini(:,1))

if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18])
stairs([0 depth(1:-1)],olddensity,'LineWidth',2)
    hold on
    stairs([0 depth_act(1:-1)],rhofirn_ini, 'LineWidth',1.5)
    plot([0, depth_weq(-1)],[c.rho_ice, c.rho_ice],'k')
    axis tight
    box on
    leg('Observations','Model initial state',\'bulk density',
        'Location','NorthOutside')
    leg boxoff
    xlabel('Depth (m)')
    ylabel('Density (kg m^{-3})','Interpreter','tex')
    view([90 90])
    title(c.station)
        print(f,sprintf('#s/Initial_rho.tif',c.OutputFolder),'-dtiff')



## ================= Initial temperature profile =========================
# same as density, we give initial temperature to the subsurface according
# to a csv file located in the Input folder
#
# See script "InitialTemperatureProfileGenerator.m" in lib folder for creation of
# the csv file.
   
# if there is thermistor record for the first time step, then reads 
# initial subsurface conditions from AWS data
oldtemp = []
time_dt = datenum(time,1,1)
clearvars out

    #reads initial subsurface conditions from file
if c.retmip == 0
    switch c.station
        case 'DYE-2_HQ'
        filename = './Input/Initial state/temperature/RetMIP_temperature_Dye-2 HQ.csv'
        case {'DYE-2','DYE-2_long'}
        filename = './Input/Initial state/temperature/RetMIP_temperature_Dye-2.csv'
        case {'Summit', 'NASA-E', 'NGRIP'}
        filename = './Input/Initial state/temperature/RetMIP_temperature_Summit.csv'
        case {'NASA-U', 'KAN-U'}
        filename = ['./Input/Initial state/temperature/RetMIP_temperature_', c.station,'.csv']
        case 'GITS'
            filename = './Input/Initial state/temperature/GITS_temperature.csv'  
        case 'TUNU-N'
            filename = './Input/Initial state/temperature/GITS_temperature.csv'  
        case 'NUK_K'
            filename = './Input/Initial state/temperature/InitialTemperatureProfile_NUK_K.csv'
        case 'Miege'
            filename = './Input/Initial state/temperature/RetMIP_temperature_Firn Aquifer.csv'

    
        otherwise
            disp('Using thermistor data')
            for i = 1:24*14 #we look at the week following the installation
                if sum(~isnan(T_obs(i,:)))>1
                    if sum(~isnan(T_obs(i,:)))<3
                        continue
                    
                    if sum(~isnan(depth_thermistor(i,~isnan(T_obs(i,:)))'))<3
                        continue
                    
                    depth= depth_thermistor(i,~isnan(T_obs(i,:)))' #in m
                    depth = depth(depth~=0)
                    [depth,ind_sorted] = sort(depth)
                    out(:,1) = depth
                    oldtemp = T_obs(i,~isnan(T_obs(i,:)))'
                    oldtemp = oldtemp(depth~=0)
                    oldtemp = oldtemp(ind_sorted)
                    out(:,2) = oldtemp
                    date_Tstring = datestr(time_dt(i))
                    break
                
            
            filename = []
        
        if isempty(oldtemp)
            disp('no thermistor data availble')
            filename = './Input/Initial state/temperature/InitialTemperatureProfile.csv'
        
        
    

else
    filename = ['.\RetMIP\Input files\temperature\RetMIP_temperature_', c.station,'.csv']


    if ~isempty(filename)
        disp('File used for firn temperature initialization:')
        disp(filename)
                delimiter = ''
        startRow = 2
        formatSpec = '#f#f#[^\n\r]'

        fileID = fopen(filename,'r')
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,\
            'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false)
        fclose(fileID)

        out = [dataArray{1:-1}]
    

# date_Tstring = 'manualy chosen profile'
depth = out(:,1) #old depth scale in m
oldtemp = out(:,2) #temperature on the old real depth scale


oldtemp(depth<0 ) =[]
depth(depth<0 ) =[]
  depth(isnan(oldtemp))=[]
  oldtemp(isnan(oldtemp))=[]
  
# Prepraring observation's depth scale
    if depth(1) ~= 0
        #if the input density profile does not contain surface temperature,
        #then we use Tsurf_in to initiate the temperature profile
        depth = [0 depth]
        oldtemp = [Tsurf_ini - c.T_0 oldtemp]
    
    if depth_act()>depth()
        depth = [depth depth_act()]
        oldtemp = [oldtemp oldtemp()]
    
        
    [newtemp] = ConvertToGivepthScale(depth, oldtemp, depth_act,'linear')

    newtemp(isnan(newtemp))=[]
    newtemp = newtemp + c.T_0 #going back to K

    # giving the observed temperature (on appropriate scale) as initial value
    # for the subsurface temperature profile
    # There might be a shift to apply deping on whether the first value in
    # subsurface column represent the temp at depth 0 (=Tsurf) or at depth 
    # c.cdel(1). To be asked.
    T_ice(1:length(newtemp),1,j) = newtemp

    # EXTRAPOLATING UNTIL DEEP FIRN TEMPERATURE
    if length(newtemp)<length(c.cdel)
        d1 = length(newtemp)
        T1 = newtemp(d1)
        d2 = c.jpgrnd
        T2 = c.Tdeep_AWS + c.T_0

        ind = length(newtemp)+1:(c.jpgrnd-1)
        T_ice(ind,1,j) = \
            (T1-T2)/(d2-d1)^2*(d2-ind).^2 + T2
        T_ice(c.jpgrnd,1,j) = c.Tdeep_AWS+ c.T_0
    
    # removing non-freezing temperatures (just in case)
    subsurfmelt = find(T_ice(:,1,j) > c.T_0)
    if sum(subsurfmelt )>0
        T_ice(subsurfmelt,1,j) = c.T_0
    
    
if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18])
    scatter(depth,oldtemp,'o','fill')
    hold on
    stairs(depth_act, T_ice(:,1,j)-c.T_0, 'LineWidth',2)
    scatter(depth(),oldtemp(),'s','LineWidth',6)
    leg('Observations','Model initial state',\
        'Prescribed deep \newlinefirn temperature',\
        'Interpreter','tex',    'Location','NorthOutside')

    ylabel('Temperature (^oC)','Interpreter','tex')
    leg boxoff
    axis tight
    box on
    view([90 90])  
    title(c.station)  
    print(f,sprintf('#s/Initial_temp',c.OutputFolder),'-dpng')


## ========== Initial grain size ===================================

switch c.station
    case 'Miege'
        out = dlmread('./Input/Initial state/InitialGrainSizeProfile_Miege.csv')
    otherwise
        out = dlmread('./Input/Initial state/InitialGrainSizeProfile.csv')


depth = out(:,1) #old depth scale in m
old_graind = out(:,2) #density in kg/m^3 on the old real depth scale

[new_graind] = ConvertToGivepthScale(depth, old_graind, depth_act,'intensive')

graind_out = 2*ones(size(c.cdel))
graind_out(1:find(~isnan(new_graind),1,'last'),1,j) = new_graind(1:find(~isnan(new_graind),1,'last'))
   
if sum(isnan(graind_out))>0
        error('NaN in initial temperature values')

if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18])
    scatter(depth,old_graind,'o')
    hold on
    stairs([0 depth_act],[graind_out graind_out()], 'LineWidth',2)
    scatter([depth_act],[graind_out], 'LineWidth',2)
    leg('Original data','Model initial state','Location','NorthOutside')

    ylabel('Grain size (mm)')
    leg boxoff
    axis tight
    xlim([0 max(depth_act)])
    box on
    view([90 90])
        print(f,sprintf('#s/Initial_graind',c.OutputFolder),'-dpng')



## ========== Initial water content ================================
switch c.station
    case 'Miege'
        [~, ind_12] = min(abs(depth_act-12))
        [~, ind_37] = min(abs(depth_act-37))
#         thick_aq = depth_act(ind_25)- depth_act(ind_15)

        pore_volume = c.cdel *1000 .* (1./rhofirn_ini - 1/c.rho_ice)
        pore_volume(snowc_ini<c.smallno) = 0
        
        saturation = 1 #2000 / (sum(pore_volume(ind_12:ind_37))*c.rho_water)
        slwc_ini(ind_12:ind_37,:) = saturation * pore_volume(ind_12:ind_37)
        fprintf('Intital total lwc of #0.2f mm\n\n', \
            1000*sum(slwc_ini(ind_12:ind_37,:)))
        
    case 'FA'
        filename = '.\RetMIP\Input files\lwc\RetMIP_lwc_FA.csv'
        delimiter = ''
        startRow = 2
        formatSpec = '#f#f#[^\n\r]'
        fileID = fopen(filename,'r')
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false)
        fclose(fileID)
        M = [dataArray{1:-1}]
        clearvars filename delimiter startRow formatSpec fileID dataArray ans
        
        depth = M(:,1)
        old_lwc = M(:,2)        
        [slwc_ini] = ConvertToGivepthScale(depth, old_lwc, depth_act,'intensive')
    otherwise
        slwc_ini = 0*depth_act
        old_lwc = 0*depth


if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18])
    stairs(depth_act,slwc_ini.*1000,'r', 'LineWidth',2)
    if c.retmip == 1
        hold on
        scatter(depth,old_lwc.*1000, 'LineWidth',2)
    
    
    hold on
    leg('Original data','Model initial state','Location','NorthOutside')
    ylabel('Initial water content (mm)')
    leg boxoff
    axis tight
    box on
#     title(sprintf('Initial grain size profile\n on #s',datestr(time_dt(1))))
    view([90 90])
    print(f,sprintf('#s/Initial_slwc.tif',c.OutputFolder),'-dtiff')


## ========== add other initialization here ========================



def [RunName, c] = OutputName(c)
if c.ConductionModel == 1
    RunName = sprintf('#s_#i_ConductionOnly',\
        tag, c.year)
else
    RunName = c.station
    RunName = [RunName, sprintf('_#i',c.year)]
    if c.calc_CLliq == 1
        text_Si = 'CL'
    else
        text_Si = sprintf('#0.2f',c.liqmax)
    
    RunName = [RunName, '_IWC_', text_Si]
    RunName = [RunName, sprintf('_#i_layers',c.jpgrnd-1)]
    
    c.OutputFolder = sprintf('#s/#s',c.OutputRoot,RunName)
    [~,~,id] =  mkdir(c.OutputFolder)
    count = 1
    while ~isempty(strfind(id,'DirectoryExists'))
       count =count+1

       c.OutputFolder = sprintf('./Output/#s_#i',RunName,count)
       [~,~,id] =  mkdir(c.OutputFolder)
    
    if count>1
           RunName = sprintf('#s_#i',RunName,count)
    


def [time, year, day, hour, pres,\
    T1, T2, z_T1, z_T2, o_T1,o_T2, \
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, \
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,\
    SRin,SRout, LRin, LRout, T_ice_obs, \
    depth_thermistor, Surface_Height, Tsurf_obs] = \
    RenamingVariables(data_out,c)

    # time
    year = data_out.Year
    hour = data_out.HourOfDayUTC
    day = data_out.DayOfYear
    # leap years    
    time = year + (day + hour/24)/365
    leapyear = find(year/4 == floor(year/4))
    if sum(leapyear) >0
        time(leapyear) = year(leapyear)+(day(leapyear)+hour(leapyear)/24.)/366
    
    
    # temperature, humidity and wind speed
    if sum(strcmp(data_out.Properties.VariableNames,'AirTemperatureC'))
        disp('Only one level was detected on the weather station.')
        T2 = data_out.AirTemperatureC
        T1 = NaN(size(T2))
        RH2 = data_out.RelativeHumidity
        RH1 = NaN(size(T2))
        WS2 = data_out.WindSpeedms
        WS1 = NaN(size(T1))

        o_T1 = NaN(size(T1))
        o_RH1 = NaN(size(T1))
        o_WS1 = NaN(size(T1))
        z_T1 = NaN(size(T1))
        z_RH1 = NaN(size(T1))
        z_WS1 = NaN(size(T1))

        # assigning height
        if sum(strcmp('HeightWindSpeedm',data_out.Properties.VariableNames))
            z_WS2 = data_out.HeightWindSpeedm
            z_T2 = data_out.HeightTemperaturem
            z_RH2 = data_out.HeightHumiditym
        else
            disp('Assigning measurement height from HeightSensorBoomm_raw field.')
            z_WS2 = data_out.HeightSensorBoomm_raw + 0.4
            z_T2 = data_out.HeightSensorBoomm_raw - 0.12
            z_RH2 = data_out.HeightSensorBoomm_raw - 0.12
        

        # assigning origin
        if sum(strcmp('WindSpeed1ms_Origin',data_out.Properties.VariableNames))
            o_WS2 = data_out.WindSpeed1ms_Origin
            o_T2 = data_out.AirTemperature1C_Origin
            o_RH2 = data_out.RelativeHumidity1_Origin
        else
            disp('No WindSpeed1ms_Origin field specified')
            o_T2 = zeros(size(T1))
            o_RH2 = zeros(size(T1))
            o_WS2 = zeros(size(T1))
        
        
    elseif sum(strcmp(data_out.Properties.VariableNames,'AirTemperature1C'))
        disp('Two levels detected on the weather station')
        T1 = data_out.AirTemperature1C
        T2 = data_out.AirTemperature2C
        RH1 = data_out.RelativeHumidity1
        RH2 = data_out.RelativeHumidity2
        WS1 = data_out.WindSpeed1ms
        WS2 = data_out.WindSpeed2ms
        
        o_T1 = data_out.AirTemperature1C_Origin
        o_T2 = data_out.AirTemperature2C_Origin
        o_RH1 = data_out.RelativeHumidity1_Origin
        o_RH2 = data_out.RelativeHumidity2_Origin
        o_WS1 = data_out.WindSpeed1ms_Origin
        o_WS2 = data_out.WindSpeed2ms_Origin

        z_T1 = data_out.HeightTemperature1m
        z_T2 = data_out.HeightTemperature2m
        z_RH1 = data_out.HeightHumidity1m
        z_RH2 = data_out.HeightHumidity2m
        z_WS1 = data_out.HeightWindSpeed1m
        z_WS2 = data_out.HeightWindSpeed2m
    else
        error('Cannot recognize temperature field in weather data file.')
    
    T1 = T1 + c.T_0       # Temperature in degrees Kelvin
    T2 = T2 + c.T_0       # Temperature in degrees Kelvin

    # radiation
    LRin = data_out.LongwaveRadiationDownWm2
    LRout = data_out.LongwaveRadiationUpWm2

    if sum(strcmp(data_out.Properties.VariableNames,'ShortwaveRadiationDown_CorWm2'))
        disp('Using ShortwaveRadiation***_CorWm2 field.')
        SRin = data_out.ShortwaveRadiationDown_CorWm2
        SRout = data_out.ShortwaveRadiationUp_CorWm2
    else
        SRin = data_out.ShortwaveRadiationDownWm2
        SRout = data_out.ShortwaveRadiationUpWm2
    

    # other variables
    pres = data_out.AirPressurehPa
    Surface_Height = data_out.SurfaceHeightm   
    Tsurf_obs = min(c.T_0, ((LRout-(1-c.em)*LRin)/c.em/c.sigma).^0.25)
    Tsurf_obs(or(isnan(LRout),isnan(LRin))) = NaN
    
    ind = strfind(data_out.Properties.VariableNames,'IceTemperature')
    ind = find(~cellfun('isempty', ind))
    ind2 = strfind(data_out.Properties.VariableNames,'DepthThermistor')
    ind2 = find(~cellfun('isempty', ind2))
    num_therm = length(ind)
    T_ice_obs = NaN(length(T1),num_therm)
    depth_thermistor = NaN(length(T1),num_therm)

    if ~isempty(ind2)
        for i = 1:length(ind)
            T_ice_obs(:,i) = data_out.(data_out.Properties.VariableNames{ind(i)})
            depth_thermistor(:,i) = data_out.(data_out.Properties.VariableNames{ind2(i)})
        
    
        

    ind =  (LRout>316)
    if sum(ind)>0
        if c.verbose == 1
        fprintf('Warning: Observed surface temperature higher than 0degC\n')
        
    #     before = LRout(ind)
    # #     LRout(ind) = LRout(ind) - (20/15 * (T(ind)-c.T_0))
    #     figure
    #     scatter(T(ind),before,'xr')
    # #     hold on
    # #     scatter(T(ind),LRout(ind),'ob')
    # #     leg('before correction', 'after correction')
    #     xlabel('Temperature (deg C)')
    #     ylabel('Outgoing long-wave radiation (W/m^2)')
    #     title('Observed LRout > black body at 0degC')
    


def [Tsurf_reset, T_ice_reset] = \
    ResetTemp(depth_thermistor, LRin, LRout, T_obs, rho, T_ice, time, k, c)

# ========== resets surface temperature ================
if ~isnan(LRout(k)) && ~isnan(LRin(k))
    Tsurf_reset = ((LRout(k) - (1-c.em)*LRin(k)) /(c.em*c.sigma))^(1/4)
else
    Tsurf_reset = NaN


# ================= resets temperature profile =========================

    # if there is thermistor record for the first time step, then reads 
    # initial subsurface conditions from AWS data
    T_ice_reset = NaN(c.jpgrnd,1)
    if sum(~isnan(T_obs(k,:)))>1
        depth = depth_thermistor(k,(depth_thermistor(k,:)~=0))'
        oldtemp = T_obs(k,(depth_thermistor(k,:)~=0))'
    
    
    # calculates the new depth scale in mweq
    depth_weq = cumsum(c.cdel)

    # calculates the new depth scale in real m
#     depth_act = depth_weq .*c.rho_water ./rho(:,k)
    depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k))

    # Here we add an initial surface temperature
    depth_act = [0 depth_act]
    depth_weq = [0 depth_weq]
    if depth(1) ~= 0
        #if the input density profile does not contain surface temperature,
        #then we use Tsurf_in to initiate the temperature profile
        depth = [0 depth]
        oldtemp = [Tsurf_reset - c.T_0 oldtemp]
    
    # the old scale is converted from m to mweq by interpolation
    oldscale_weq = interp1(depth_act,depth_weq,depth)

    # we limit the new scale to the values that are given within the oldscale
    newscale_weq = depth_weq(depth_weq<=oldscale_weq())

    # the temperature is interpolated on each stamp of the new scale
    newtemp = interp1(oldscale_weq,oldtemp,newscale_weq)
    newtemp = newtemp + c.T_0 #going back to K

    # giving the observed temperature (on appropriate scale) as initial value
    # for the subsurface temperature profile
    # There might be a shift to apply deping on whether the first value in
    # subsurface column represent the temp at depth 0 (=Tsurf) or at depth 
    # c.cdel(1). To be asked.
    T_ice_reset(1:length(newtemp)) = newtemp

    # fills the rest of the temperature profile withan interpolation that
    # s with Tdeep at the bottom of the profile

    if length(newtemp)<length(c.cdel)
        d1 = length(newtemp)
        T1 = newtemp(d1)
        d2 = c.jpgrnd
        T2 = c.Tdeep_AWS + c.T_0
        ind = length(newtemp)+1:(c.jpgrnd-1)
        T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) = \
            (T1-T2)/(d2-d1)^2*(d2-ind).^2 + T2
        T_ice_reset(c.jpgrnd) = c.Tdeep_AWS+ c.T_0

#         T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) = NaN
    
    
#     figure
#     scatter(oldscale_weq,oldtemp,'o')
#     hold on
#     stairs(depth_weq(1:-1),T_ice(:,k,1)-c.T_0, 'LineWidth',2)
#     stairs(depth_weq(1:-1),T_ice_reset-c.T_0, 'LineWidth',2)
#     leg('data','before reset','after reset','Location','South')
#     xlabel('Depth (m weq)')
#     ylabel('Temperature (deg C)')
#     title(sprintf('2/ #s',datestr(datenum(time(k),0,0))))
#     view([90 90])
    
# removing non-freezing temperatures (just in case)
subsurfmelt = find(T_ice_reset(:) > c.T_0)
if sum(subsurfmelt )>0
    T_ice_reset(subsurfmelt) = c.T_0

## Comes from outside the def just after it
#                     figure
                    depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k))
                    depth_act = [0 depth_act]

#                     scatter(depth_thermistor(k,depth_thermistor(k,:)~=0),\
#                         T_ice_obs(k,depth_thermistor(k,:) ~= 0), 'o')
#                     hold on
#                     stairs(depth_act(1:-1),T_ice(:,k,j)-c.T_0)
                    
                T_ice(~isnan(T_reset),k,j) = T_reset(~isnan(T_reset))
#                     stairs(depth_act(1:-1),T_ice(:,k,j)-c.T_0)
#                     leg('data','before reset','after reset','Location','South')
#                     xlabel('Depth (m)')
#                     ylabel('Temperature (deg C)')
#                     title(sprintf('#s',datestr(datenum(time(k),0,0))))
#                     view([90 90])
   
                    [zso_capa, zso_cond] = ice_heats (c)
                    [grndc, grndd, ~, ~]\
                        = update_tempdiff_params (rho(:,k), Tdeep(j)                    \
                        , snowc, snic, T_ice(:,k,j), zso_cond, zso_capa, c)



def RH_wrt_w = RHice2water(RH_wrt_i, T, pres)
# def converting humidity with regards to ice into humidity with
# regards to water using GOFF-GRATCH (1945) formula.
#   
#   RH_wrt_i is an single value or an array of relative humidity with
#   regards to ice given in percent
#
#   T is the corresponding air temperature in degree Celsius
#
#   pres is the corresponding air pressure in hPa (not used in the current
#   form of the def)

T_0 = 273.15
T_100 = 373.15
# eps is the ratio of the molecular weights of water and dry air
eps = 0.62198 

# Rv = 461.5
# es_wtr = eps * exp( 1/Rv * ((L + T_0 * beta)*(1/T_0 - 1/T) - beta* log(T./T_0)))

# GOFF-GRATCH 1945 equation
es_wtr = 10.^( \
    -7.90298*(T_100./T - 1) + 5.02808 * log10(T_100./T) \   # saturation vapour pressure above 0 C (hPa)
    - 1.3816E-7 * (10.^(11.344*(1.-T./T_100))-1.) \
    + 8.1328E-3*(10.^(-3.49149*(T_100./T-1)) -1.) + log10(1013.246) )

# WMO 2012 equation (based on Goff 1957)
# es_wtr = 10.^( \
#     10.79574*(1 - T_0./T) + 5.028 * log10(T / T_0) \   # saturation vapour pressure above 0 C (hPa)
#     + 1.50475E-4 * (1 - 10.^(-8.2969 * (T./T_0 - 1))) \
#     + 0.42873E-3*(10.^(4.76955*(1 - T_0./T)) -1.) +  0.78614 + 2.0 )

es_ice = 10.^( \
    -9.09718 * (T_0 ./ T - 1.) - 3.56654 * log10(T_0 ./ T) + \
    0.876793 * (1. - T ./ T_0) + log10(6.1071)  )   # saturation vapour pressure below 0 C (hPa)

# es_ice = 10.^( \
#     -9.09685 * (T_0 ./ T - 1.) - 3.56654 * log10(T_0 ./ T) + \
#     0.87682 * (1. - T ./ T_0) + 0.78614 + 2.0  ) 

# q_sat_wtr = eps * es_wtr./(pres-(1-eps)*es_wtr)    # specific humidity at saturation over water
# q_sat = eps * es_ice./(pres-(1-eps)*es_ice)        # specific humidity at saturation over ice
freezing=ones(size(T))
freezing(T>=T_0) = 0
freezing = freezing==1

RH_wrt_w = RH_wrt_i
RH_wrt_w(freezing) = RH_wrt_i(freezing) .*es_ice(freezing) ./ es_wtr (freezing)               # specific humidity in kg/kg
RH_wrt_w(RH_wrt_w<0)=0
RH_wrt_w(RH_wrt_w>100) = 100



def [RH_wrt_i, RH_wrt_w] = spechum2relhum(T, pres, q, c)

# spechum2relhum
#
# Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
# translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
#==========================================================================

    # SPECIFIC HUMIDITY & SATURATION -----------------------------------------------------------------------
    es_wtr = 10.^(-7.90298*(c.T_100./T-1) + 5.02808 * log10(c.T_100./T) \   # saturation vapour pressure above 0 C (hPa)
        - 1.3816E-7 * (10.^(11.344*(1.-T./c.T_100))-1.) \
        + 8.1328E-3*(10.^(-3.49149*(c.T_100./T-1)) -1.) + log10(c.es_100))

    es_ice = 10.^(-9.09718 * (c.T_0 ./ T - 1.) - 3.56654 * log10(c.T_0 ./ T) + \
        0.876793 * (1. - T ./ c.T_0) + log10(c.es_0))   # saturation vapour pressure below 0 C (hPa)
    
    q_sat = c.es * es_wtr./(pres-(1-c.es)*es_wtr)    # specific humidity at saturation (incorrect below melting point)

    freezing = find(T < c.T_0)            # replacing saturation specific humidity values below melting point
    RH_wrt_w = q ./ q_sat*100
    if sum(freezing) > 0
        q_sat(freezing) = c.es * es_ice(freezing)./(pres(freezing)-(1-c.es)*es_ice(freezing))
    

    RH_wrt_i = q ./ q_sat*100
    supersaturated = find(RH_wrt_i > 100)       # replacing values of supersaturation by saturation

    if sum(supersaturated) > 0
        RH_wrt_i(supersaturated) = 100
    


def [c] = ImportConst(param)
# ImportConst: Reads physical, site-depant, simulation-depant and
# user-defined parameters from a set of csv files located in the ./Input
# folder. It stores all of them in the c structure that is then passed to
# all defs. Structures were found to be the fastest way of
# communicating values from one def to another.
#
# Author: Baptiste Vandecrux (bava@byg.dtu.dk)
# ========================================================================

## Import constants for the transect-mode run ----------------------------
# originally the surface energy balance was designed to work on transects.
# This defnality is not working anymore but might be implemented again
# later on.

filename = 'const_transect.csv'
delimiter = ''
startRow = 2
formatSpec = '#s#f#s#[^\n\r]'
fileID = fopen(filename,'r')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false)
fclose(fileID)
Parameter = dataArray{:, 1}
Value = dataArray{:, 2}
clearvars filename delimiter startRow formatSpec fileID dataArray ans

for i=1:length(Parameter)
    eval(sprintf('c.#s=#f',Parameter{i},Value(i)))


clear Parameter Value
c.rows = 176038

## Import constants regarding the station --------------------------------
filename = 'StationInfo.csv'
delimiter = ''
startRow = 2
formatSpec = '#s#f#f#f#f#f#f#f#f#[^\n\r]'
fileID = fopen(filename,'r')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false)
fclose(fileID)
StationInfo = table(dataArray{1:-1}, 'VariableNames', {'stationname','latitude','longitude','elevationm','deepfirntemperaturedegC','slopedeg','meanaccumulationm_weq','InitialheightTm','InitialheightWSm'})
clearvars filename delimiter startRow formatSpec fileID dataArray ans

i_station = find(strcmp(param.station, StationInfo.stationname))
c.lat = StationInfo.latitude(i_station)
c.lon = StationInfo.longitude(i_station)
c.elev_AWS = StationInfo.elevationm(i_station)
c.ElevGrad = StationInfo.slopedeg(i_station)*pi/180 # converting slope in degrees to elevation gradient (dz/dx) in m/m positive downward
c.H_T = StationInfo.InitialheightTm(i_station)
c.H_WS = StationInfo.InitialheightWSm(i_station)
c.Tdeep_AWS = StationInfo.deepfirntemperaturedegC(i_station)
c.accum_AWS = StationInfo.meanaccumulationm_weq(i_station)

## Import simulation constant -----------------------------------------------------------------------

filename = 'const_sim.csv'
delimiter = ''
startRow = 1
formatSpec = '#s#s#s#[^\n\r]'
fileID = fopen(filename,'r')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false)
fclose(fileID)
Parameter = dataArray{:, 1}
Value = dataArray{:, 2}
clearvars filename delimiter startRow formatSpec fileID dataArray ans

for i=1:length(Parameter)
    if isempty(Parameter{i})
        continue
    
    eval(sprintf('c.#s=#s',Parameter{i},Value{i}))

clear Parameter Value

## Import physical constants -----------------------------------------------------------------------
filename = 'const_phy.csv'
delimiter = {',',''}
formatSpec = '#s#s#s#s#s#[^\n\r]'
fileID = fopen(filename,'r')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false)
fclose(fileID)
raw = repmat({''},length(dataArray{1}),length(dataArray)-1)
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col}

numericData = NaN(size(dataArray{1},1),size(dataArray,2))

for col=[2,3,4]
    # Converts strings in the input cell array to numbers. Replaced non-numeric
    # strings with NaN.
    rawData = dataArray{col}
    for row=1:size(rawData, 1)
        # Create a regular expression to detect and remove non-numeric prefixes and
        # suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)'                result = regexp(rawData{row}, regexstr, 'names')
        
        try
            numbers = result.numbers
            
            # Detected commas in non-thousand locations.
            invalidThousandsSeparator = false
            if any(numbers==',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$'
                if isempty(regexp(thousandsRegExp, ',', 'once'))
                    numbers = NaN
                    invalidThousandsSeparator = true
                
            
            # Convert numeric strings to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(strrep(numbers, ',', ''), '#f')
                numericData(row, col) = numbers{1}
                raw{row, col} = numbers{1}
            
        catch
        
        
    


rawNumericColumns = raw(:, [2,3,4])
rawCellColumns = raw(:, [1,5])
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns) # Find non-numeric.lls
rawNumericColumns(R) = {NaN} # Replace non-numeric.lls
Parameter = rawCellColumns(:, 1)
Value1 = cell2mat(rawNumericColumns(:, 1))
Value2 = cell2mat(rawNumericColumns(:, 2))
Value3 = cell2mat(rawNumericColumns(:, 3))
# clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R

for i=1:length(Parameter)
    switch Parameter{i}
        case {'ch1','ch2','ch3','cq1','cq2','cq3'}
            eval(sprintf('c.#s=[#f #f #f]',Parameter{i},Value1(i),Value2(i),Value3(i)))
        otherwise
            eval(sprintf('c.#s = #e',Parameter{i},Value1(i)))
    

# clear Parameter Value1 Value2 Value3

## Import constants used by the subsurface model -----------------
filename = 'const_subsurf.csv'

delimiter = ''
startRow = 1
formatSpec = '#s#s#s#[^\n\r]'
fileID = fopen(filename,'r')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false)
fclose(fileID)
Parameter = dataArray{:, 1}
Value = dataArray{:, 2}
clearvars filename delimiter startRow formatSpec fileID dataArray ans

for i=1:length(Parameter)
    if ~isempty(Parameter{i})
#         disp(sprintf('c.#s=#s',strtrim(Parameter{i}),Value{i}))
        eval(sprintf('c.#s=#s',strtrim(Parameter{i}),Value{i}))
    

clear Parameter Value

## Other global parameters

# Update B2017 now the thickness of a new layer is set to the fifth of the
# mean annual precipiation amount.
c.lim_new_lay = c.accum_AWS/c.new_lay_frac


varname1 = fieldnames(c)
varname2 = fieldnames(param)
if param.verbose == 1
    fprintf('\nOverwriting default value for:\n')

for i = 1:length(varname1)
    for j =1:length(varname2)
        if strcmp(varname1(i),varname2(j))
            c.(char(varname1(i))) = param.(char(varname2(j)))
        
    

if c.verbose == 1
for i = 1:length(varname1)
    for j =1:length(varname2)
        if strcmp(varname1(i),varname2(j))           
            fprintf('#s #s #0.2e\n',char(varname1(i)),\
                repmat(' ',1,20-length(char(varname1(i)))),\
                param.(char(varname2(j))))
        
    


c.InputAWSFile = param.InputAWSFile
c.station = param.station

if isfield(param,'cdel') #if cdel has been defined in param
    c.cdel = param.cdel
    c.z_ice_max     = length(c.cdel)-1   # number of sub-surface levels
    c.jpgrnd = c.z_ice_max+1
else
    c.z_ice_max     = c.z_max/c.dz_ice   # number of sub-surface levels
    c.jpgrnd = c.z_ice_max+1
#         c.cdel = ones(c.jpgrnd,1)*c.dz_ice
    thick_cumul =  (1:c.jpgrnd).^4/c.jpgrnd^4 *c.z_max
    c.cdel = [thick_cumul(1) thick_cumul(2:) - thick_cumul(1:-1)]'
    tmp = max(0,c.lim_new_lay - c.cdel)
    c.cdel = c.cdel + tmp - flipud(tmp)


c.rh2oice = c.rho_water/c.rho_ice
c.cmid = zeros(c.jpgrnd,1)
c.rcdel = zeros(c.jpgrnd,1)

c.cdelsum = 0
for jk = 1:c.jpgrnd
    c.cdelsum = c.cdelsum + c.cdel(jk)
    c.cmid(jk) = c.cdelsum - ( c.cdel(jk) / 2 )
    c.rcdel(jk) = 1/c.cdel(jk)

c.cdelV      =zeros(c.jpgrnd,1)
c.zdtime = c.delta_time

# Determine local runoff time-scale  (Zuo and Oerlemans 1996). Parameters
# are set as in Lefebre et al (JGR, 2003) = MAR value (Fettweis pers comm)
c.t_runoff = (c.cro_1 + c.cro_2 * exp(- c.cro_3 * c.ElevGrad))*c.t_runoff_fact





def [var_new] = ConvertToGivepthScale(depth_old, var_old, depth_new,opt)
# Interpolates the depth profile of a given variable old_var 
# (temperature grain size\.) according to a given scale new_depth in real m

    transpose_at_the_ = 0

    if isrow(depth_old) ~= isrow(depth_new)
        error('Old and new scale should be both rows or both columns.')
    elseif ~isrow(depth_old)
        transpose_at_the_ = 1
       depth_old=depth_old'
       var_old=var_old'
       depth_new = depth_new'
    
    var_new = NaN(size(depth_new))

    switch opt
        case 'linear'
        # the variable is interpolated on each stamp of the new scale
        var_new = interp1(depth_old,var_old,depth_new,'linear','extrap')
        
        case 'intensive'
            # intensive variables do not dep on the size of the system
            # f.e. if you know the density of a core section and you cut it
            # in two, the two sub-sections can be assigned the same density
            # as the original section. However we need to take into account
            # into the re-sampling the case when a new section is composed
            # of two sections in the old scale. Then the density of that
            # new section is the thickness-weighted average of the two
            # original sections.
            # example: depth_old = 1:4 density_old = 100:100:400
            #           depth_new = [0.1 0.2 0.6 1.2 3.5]  
            # density_new = [100.0000  100.0000  100.0000  133.3333
            # 286.9565]

            if depth_new()>depth_old()
                depth_old = [depth_old, depth_new()]
                var_old = [var_old, var_old()]
            
            left_neighbour_in_new = depth_new
            left_neighbour_in_new(1) = 0
            left_neighbour_in_new(2:) = depth_new(1:-1)

            left_neighbour_in_old =  interp1([0 depth_old],[0 depth_old],depth_new,'previous')
            
            ind_type_1 = left_neighbour_in_new >= left_neighbour_in_old
            ind_type_2 = find(left_neighbour_in_new < left_neighbour_in_old)
            
            var_new(ind_type_1) = interp1([0 depth_old],\
                [var_old(1) var_old],\
                depth_new(ind_type_1),'next')
            
            depth_merged = [depth_old depth_new]
            var_merged = [var_old var_new]
            
            [depth_merged, ind_sorted] = sort(depth_merged)
            var_merged = var_merged(ind_sorted)
            var_merged(isnan(var_merged)) = interp1([0  depth_merged(~isnan(var_merged))],\
                [var_merged(1) var_merged(~isnan(var_merged))],\
                depth_merged(isnan(var_merged)),'next')

            thick_merged = depth_merged
            thick_merged(2:) = depth_merged(2:) - depth_merged(1:-1)
            
            for i = ind_type_2
                i_in_merged = discretize(depth_merged, \
                    [left_neighbour_in_new(i) depth_new(i)])
                i_in_merged(isnan(i_in_merged))= 0
                if i~=1
                    i_in_merged(find(i_in_merged,1,'first')) = 0 # transforming the first one into 0
                
                var_new(i) = sum(var_merged(i_in_merged==1).*thick_merged(i_in_merged==1)) \
                    ./sum(thick_merged(i_in_merged==1))
            
            
        case 'extensive'
            # extensive values dep on the size of the system
            # for example when downscaling the liquid water content of one
            # cell into two equally sized cells then the lwc in the new
            # cells are half of the lwc of the original ones
            # example: depth_old = 1:4 lwc_old = [1 0 0 1]
            #           depth_new = [0.1 0.2 0.6 1.2 3.5]  
            
            if depth_new()>depth_old()
                thick_last_old = depth_old()-depth_old(-1)
                depth_old = [depth_old, depth_new()]
                thick_last_old_new = depth_old()-depth_old(-1)
                var_old = [var_old, var_old()/thick_last_old*thick_last_old_new]
            
            
            depth_merged = sort([depth_old depth_new])
                        
            thick_merged = depth_merged
            thick_merged(2:) = depth_merged(2:) - depth_merged(1:-1)
            
            thick_old = depth_old
            thick_old(2:) = depth_old(2:) - depth_old(1:-1)

            ind_bin = discretize(depth_merged,[0 depth_old],'IncludedEdge','right')

            if sum(isnan(ind_bin))>1
                error('Some depths asked in depth_new not covered by depth_old')
            
            var_merged = var_old(ind_bin).* thick_merged ./ thick_old(ind_bin)
            
            ind_bin_new = discretize(depth_merged,[0 depth_new],'IncludedEdge','right')
            var_merged(isnan(ind_bin_new)) =  []
            ind_bin_new(isnan(ind_bin_new)) =  []
            
            var_new = accumarray(ind_bin_new', var_merged')'
    
    
    if transpose_at_the_==1
       var_new=var_new'
    


