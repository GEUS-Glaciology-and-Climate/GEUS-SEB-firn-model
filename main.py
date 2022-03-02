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
snic[:, -1] = df_ini.snic
snowc[:, -1] = df_ini.snowc
dgrain[:, -1] = df_ini.grain_size_mm
tsoil[:, -1] = df_ini.temp_degC
grndc[:, -1] = tsoil[:, -1]
snowbkt[-1] = 0
i = 0

# %% processing step by step
# print(i)
# (pts, pgrndc, pgrndd, pslwc, psnic, psnowc, prhofirn, ptsoil, pdgrain, zsn, zraind, zsnmel, pTdeep, psnowbkt, c) = (ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(), slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(),  rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(), net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)

# for i in progressbar(range(0, len(time))):
# snowc[:,i], snic[:,i],  slwc[:,i], tsoil[:,i],  zrfrz[:,i], \
# rhofirn[:,i], zsupimp[:, i],  dgrain[:,i],  zrogl[i],  ts_out[i], \
# grndc[:,i], grndd[:,i],  pgrndcapc[i], pgrndhflx[i], dH_comp[i], \
# snowbkt[i], compaction[:,i] = \
#     sub.subsurface(ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
#                     slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(), 
#                     rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
#                     net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)

# depth_weq_zero = sub.numba_insert(psnowc, 0, 0)
# ptsoil = sub.tsoil_diffusion(pts, pgrndc, pgrndd, ptsoil)
# prhofirn, dH_comp, compaction = sub.densification(pslwc, psnowc, psnic, prhofirn, ptsoil, c)


# i = i+1

# %% processing all
for i in progressbar(range(0, len(time))):
    snowc[:,i], snic[:,i],  slwc[:,i], tsoil[:,i],  zrfrz[:,i], \
    rhofirn[:,i], zsupimp[:, i],  dgrain[:,i],  zrogl[i],  ts_out[i],  \
    grndc[:,i], grndd[:,i],  pgrndcapc[i], pgrndhflx[i], dH_comp[i],  \
    snowbkt[i], compaction[:,i] = \
        sub.subsurface(ts[i].copy(), grndc[:,i-1].copy(), grndd[:, i-1].copy(),
                       slwc[:, i-1].copy(), snic[:, i-1].copy(), snowc[:, i-1].copy(), 
                       rhofirn[:, i-1].copy(), tsoil[:, i-1].copy(), dgrain[:, i-1].copy(),
                       net_accum[i].copy(), 0, melt[i].copy(), c.Tdeep, snowbkt[i-1].copy(), c)
# 100% (4958 of 4958) |####################| Elapsed Time: 0:00:20 Time:  0:00:20
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
