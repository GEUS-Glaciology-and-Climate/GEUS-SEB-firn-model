# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd

df = pd.read_csv('input/weather data/CARRA_vs_AWS.csv')

print(df.groupby('var')[['ME','RMSE']].mean())
