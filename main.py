# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import numpy as np
import lib_subsurface as sub
import lib_initialization as ini
from progressbar import progressbar

c = ini.ImportConst()
c.station = 'IMAU_aws4'
c.rh2oice = c.rho_water/c.rho_ice
c.zdtime = 3600
c.ElevGrad = 0.1
NumLayer = 51
c.num_lay = NumLayer
c.z_max = 50
c.dz_ice = c.z_max/NumLayer
c.verbose = 1

c.Tdeep = 250.15
# c.lim_new_lay = c.accum_AWS/c.new_lay_frac;
c.lim_new_lay = 0.02
c.rho_fresh_snow = 315

df_aws = pd.read_csv('./Input/'+c.station+'_high-res_meteo.csv', sep=';')
df_aws.time = pd.to_datetime(df_aws.time)
df_aws = df_aws.set_index('time').resample('H').mean()

time = df_aws.index.values
ts = df_aws.Tsurf_K.values
net_accum = df_aws.acc_subl_mmweq.values/1000
melt = df_aws.melt_mmweq.values/1000

df_ini = ini.InitializationSubsurface(c)

rhofirn = np.empty((c.num_lay, len(time)))
snowc = np.empty((c.num_lay, len(time),))
snic = np.empty((c.num_lay, len(time),))
slwc = np.empty((c.num_lay, len(time),))
dgrain = np.empty((c.num_lay, len(time),))
tsoil = np.empty((c.num_lay, len(time),))
grndc = np.empty((c.num_lay, len(time),))
grndd = np.empty((c.num_lay, len(time),))
compaction = np.empty((c.num_lay, len(time),))
zrfrz = np.empty((c.num_lay, len(time),))
zsupimp = np.empty((c.num_lay, len(time),))

ts_out = np.empty((len(time),))
zrogl = np.empty((len(time),))
pgrndcapc = np.empty((len(time),))
pgrndhflx = np.empty((len(time),))
dH_comp = np.empty((len(time),))
snowbkt = np.empty((len(time)))

# first time step
rhofirn[:, -1] = df_ini.rhofirn
# rhofirn[:, -1] = np.array([313., 313.27634074, 313.63636364, 318.17654149, 316.70683053, 311., 317.39479791, 322.57618704, 343.63636364, 315.20300604, 300., 300., 309.58487339, 300., 344.71467171, 334.37708337, 331.30321981, 413.42150125, 366.52149677, 368.32330765, 379.18583305, 413.34270811, 464.55351657, 452.71846884, 462.53169247, 477.95, 493.28730491, 504.08349529, 520.20355657, 533.33505172, 544.85678291, 557.78742916, 569.57469473, 582.00788382, 590.98335659, 605.86596393, 616.68262701, 629.67858543, 644.46083704, 656.06651191, 675.6781868, 689.42798196, 700.68471921, 719.51968441, 733.01163776, 750.43080583, 761.71887162, 843.18484452, 856.1015172, 864.77253642, 868.26566916])

snic[:, -1] = df_ini.snic
snowc[:, -1] = df_ini.snowc
dgrain[:, -1] = df_ini.grain_size_mm
# dgrain[:, -1] = [0.202464598259021, 0.554759797823622, 0.594095564235787, 0.603829864237221, 0.611536223720032, 0.620056658078672, 0.701513049031701, 0.883058339142298, 0.986753792811423, 0.994636022215787, 1.006302308412469, 1.028502121070386, 1.185929391981719, 1.254225777454145, 1.838551462924819, 1.832905967507535, 1.879271351823968, 2.048582286513246, 2.074835754536944, 2.124769392864410, 2.205411941679853, 2.273737991432119, 2.384466688513645, 2.375913292230251, 2.350407249392486, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000, 2.343000000000000]
tsoil[:, -1] = df_ini.temp_degC + c.T_0
# tsoil[:, -1] = 100 * np.array([2.454632199746392, 2.452866121740167, 2.451125244549221, 2.449376288636774, 2.447595239492818, 2.445850074520300, 2.444132941309519, 2.442521044716536, 2.440287794311223, 2.437121247519070, 2.432962638879085, 2.427788951176409, 2.424298739136632, 2.420582331305712, 2.415903113015880, 2.412317335354643, 2.410361255712495, 2.408399680000000, 2.408399680000000, 2.409581295908707, 2.411630805639346, 2.413226054628521, 2.415701003436864, 2.418121253739697, 2.419386272885430, 2.419628342997080, 2.419673878962461, 2.419722993528668, 2.419776119399773, 2.419833590466084, 2.419895433067067, 2.419961947757832, 2.420033235964409, 2.420109917119185, 2.420191408194301, 2.420278427130226, 2.420370847452908, 2.420468565392353, 2.420572229496365, 2.420680723755696, 2.420795125809470, 2.420915929051683, 2.421041781694274, 2.421173839008584, 2.421311604234549, 2.421456414376028, 2.421500000000000, 2.421500000000000, 2.421500000000000, 2.421500000000000, 2.5015000000000005])
grndc[:, -1] = tsoil[:, -1]
snowbkt[-1] = 0
i = 0

# %% processing step by step
# print(i)
# (pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep, psnowbkt, c) = (ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(), slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(),  rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(), net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)

# # for i in progressbar(range(0, len(time))):
# snowc[:,i], snic[:,i],  slwc[:,i], tsoil[:,i],  zrfrz[:,i], \
# rhofirn[:,i], zsupimp[:, i],  dgrain[:,i],  zrogl[i],  ts_out[i], \
# grndc[:,i], grndd[:,i],  pgrndcapc[i], pgrndhflx[i], dH_comp[i], \
# snowbkt[i], compaction[:,i] = \
#     sub.subsurface(ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
#                    slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(), 
#                    rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
#                    net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)

# i = i+1

# %% processing all

for i in progressbar(range(0, len(time))):
    snowc[:,i], snic[:,i],  slwc[:,i], tsoil[:,i],  zrfrz[:,i], \
    rhofirn[:,i], zsupimp[:, i],  dgrain[:,i],  zrogl[i],  ts_out[i], \
    grndc[:,i], grndd[:,i],  pgrndcapc[i], pgrndhflx[i], dH_comp[i], \
    snowbkt[i], compaction[:,i] = \
        sub.subsurface(ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
                       slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(), 
                       rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
                       net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)

# %% 
thickness_act = snowc * (c.rho_water / rhofirn)+ snic * (c.rho_water / c.rho_ice)
depth_act = np.cumsum(thickness_act, 0)
density_bulk = (snowc + snic) / (snowc / rhofirn + snic / c.rho_ice)
#%% writing output
import xarray as xr
import lib_io as io
import os
c.RunName = c.station + '_' + str(c.num_lay)+'_layers'
i = 0
succeeded = 0
while succeeded == 0:
    try:
        os.mkdir('./Output/'+c.RunName)
        succeeded = 1
    except:
        if i == 0:
            c.RunName = c.RunName + '_'+str(i)
        else:
            c.RunName = c.RunName[:-len(str(i-1))]+str(i)
        i = i+1
        
io.write_2d_netcdf(snowc, 'snowc', depth_act, time, c)
io.write_2d_netcdf(snic, 'snic', depth_act, time, c)
io.write_2d_netcdf(slwc, 'slwc', depth_act, time, c)
io.write_2d_netcdf(density_bulk, 'density_bulk', depth_act, time, c)
io.write_2d_netcdf(rhofirn, 'rhofirn', depth_act, time, c)
io.write_2d_netcdf(tsoil, 'T_ice', depth_act, time, c)
# io.write_2d_netcdf(rhofirn, 'rho_firn_only', depth_act, time, RunName)
# io.write_2d_netcdf(rfrz, 'rfrz', depth_act, time, RunName)
# io.write_2d_netcdf(dgrain, 'dgrain', depth_act, time, RunName)
# io.write_2d_netcdf(compaction, 'compaction', depth_act, time, RunName)
print(c.RunName)
