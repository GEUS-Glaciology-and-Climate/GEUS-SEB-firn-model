# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Runs the model of Covi et al. (2022) (https://github.com/fcovi/1D_SEB_model,
cloned in ../Covi_1D_SEB_model) at KAN_U for 2023, without modifying its code:
the parameters of its params.py are overwritten before the model is imported.

Needs the Python 3.9 / pandas 1.5 environment the model was written for:
    conda run -n covi_seb python tools/covi_kanu_2023/run_covi.py skin
    conda run -n covi_seb python tools/covi_kanu_2023/run_covi.py observed

Tsurf modes:
    skin      surface temperature from the skin layer formulation (sim_ref)
    observed  surface temperature from measured LWout (sim_Tsobs)
"""
import os
import sys
import time

covi_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         '../../../Covi_1D_SEB_model'))
mode = sys.argv[1]

os.chdir(covi_path)
sys.path.insert(0, covi_path)
import params  # noqa: E402  (the model's own parameter module)

params.path = covi_path
params.run_name = f'KAN_U2023_{mode}'
params.met_file = f'{covi_path}/input/KAN_U_2023/KAN_U2023_metinit.txt'
params.firn_file = f'{covi_path}/input/KAN_U_2023/KAN_U2023_firninit_5cm35m.txt'
params.radpen_folder = f'{covi_path}/radpen'
params.elev = 1840          # KAN_U
params.prec_correction = 1  # snowfall already derived from the SR50 record
params.method_Tsurf = {'skin': 2, 'observed': 1}[mode]

from core import Model  # noqa: E402

start_time = time.time()
model = Model()
model.model()

os.makedirs('output', exist_ok=True)
model.surf.to_csv(f'output/surface_{params.run_name}.txt', float_format='%.6f', na_rep='nan')
model.subs.to_dataframe().to_csv(f'output/subsurface_{params.run_name}.txt',
                                 float_format='%.6f', na_rep='nan')
with open(f'output/params_{params.run_name}.txt', 'w') as f:
    for k in dir(params):
        if not k.startswith('__'):
            f.write(f'{k} = {getattr(params, k)!r}\n')
print('--- %s seconds ---' % (time.time() - start_time))
