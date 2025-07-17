# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pickle 
with open('output/spin up 3H/Summit_100_layers_3H/Summit_final.pkl', 'rb') as f:
    out = pickle.load(f)

names = [
    'snowc', 'snic', 'slwc', 'T_ice',  'rhofirn',
    'dgrain', 'Tsurf', 'grndc', 'grndd', 'snowbkt', 
]
import matplotlib.pyplot as plt
for i in range(len(out)):
    plt.figure()
    plt.plot(out[i])
    plt.ylabel(names[i])