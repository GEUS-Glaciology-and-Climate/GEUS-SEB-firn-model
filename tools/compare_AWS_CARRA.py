import __init__
import lib_io as io
import numpy as np
import pandas as pd

from lib_initialization import load_json
from lib_CARRA_initialization import load_CARRA_data_opt

# Function to compare the output data from AWS and CARRA
def compare_data():
    length = 4373 * 3

    # Data from AWS
    # DataFrame with the weather data is created
    # Read paths for weather input and output file
    parameters = load_json()
    weather_station = str(parameters['weather_station'])
    weather_data_input_path_unformatted = str(parameters['weather_data']['weather_input_path'])
    weather_data_input_path = weather_data_input_path_unformatted.format(weather_station)

    df_aws = io.load_promice(weather_data_input_path)[3858:(3858 + 8999)]
    df_aws = df_aws.set_index("time").resample("H").mean()
    df_aws = df_aws.interpolate()
 
    time_AWS = df_aws.index.values
    T_AWS = df_aws.AirTemperature1C.values + 273.15
    z_T_AWS = df_aws.HeightTemperature1m.values
    RH_AWS = df_aws.RelativeHumidity1.values
    z_RH_AWS = df_aws.HeightHumidity1m.values
    WS_AWS = df_aws.WindSpeed1ms.values
    z_WS_AWS = df_aws.HeightWindSpeed1m.values
    pres_AWS = df_aws.AirPressurehPa.values
    SRin_AWS = df_aws.ShortwaveRadiationDownWm2.values
    SRout_AWS = df_aws.ShortwaveRadiationUpWm2.values
    LRin_AWS = df_aws.LongwaveRadiationDownWm2.values
    LRout_AWS = df_aws.LongwaveRadiationUpWm2.values
    snowfall_AWS = df_aws.Snowfallmweq.values
    rainfall_AWS = df_aws.Rainfallmweq.values

    # 3hourly AWS data time
    i_AWS = 0
    third_len_AWS = int(len(time_AWS)/3)
    time_AWS3 = np.copy(time_AWS[:third_len_AWS])
    T_AWS3 = np.copy(T_AWS[:third_len_AWS])
    z_T_AWS3 =np.copy(T_AWS[:third_len_AWS])
    RH_AWS3 = np.copy(T_AWS[:third_len_AWS])
    z_RH_AWS3 = np.copy(T_AWS[:third_len_AWS])
    WS_AWS3 = np.copy(T_AWS[:third_len_AWS])
    z_WS_AWS3 = np.copy(T_AWS[:third_len_AWS])
    pres_AWS3 = np.copy(T_AWS[:third_len_AWS])
    SRin_AWS3 = np.copy(T_AWS[:third_len_AWS])
    SRout_AWS3 = np.copy(T_AWS[:third_len_AWS])
    LRin_AWS3 = np.copy(T_AWS[:third_len_AWS])
    LRout_AWS3 = np.copy(T_AWS[:third_len_AWS])
    snowfall_AWS3 = np.copy(T_AWS[:third_len_AWS])
    rainfall_AWS3 = np.copy(T_AWS[:third_len_AWS])

    for i in range(0, len(time_AWS3)):
        time_AWS3[i] = time_AWS[i_AWS]
        T_AWS3[i] = T_AWS[i_AWS]
        z_T_AWS3[i] = z_T_AWS[i_AWS]
        RH_AWS3[i] = RH_AWS[i_AWS]
        z_RH_AWS3[i] = z_RH_AWS[i_AWS] 
        WS_AWS3[i] = WS_AWS[i_AWS]
        z_WS_AWS3[i] = z_WS_AWS[i_AWS]
        pres_AWS3[i] = pres_AWS[i_AWS]
        SRin_AWS3[i] = SRin_AWS[i_AWS]
        SRout_AWS3[i] = SRout_AWS[i_AWS]
        LRin_AWS3[i] = LRin_AWS[i_AWS]
        LRout_AWS3[i] = LRout_AWS[i_AWS]
        snowfall_AWS3[i] = snowfall_AWS[i_AWS]
        rainfall_AWS3[i] = rainfall_AWS[i_AWS]
        i_AWS = i_AWS + 3

    start_date_AWS = time_AWS3[0]
    end_date_AWS = time_AWS3[-1]
    print("AWS, start date:", start_date_AWS, ", end date: ", end_date_AWS)
    start_index = 54327
    end_index = 54327 + (len(time_AWS))
 
    weather_station = 'KAN_U'
    #weather_df = load_CARRA_data(weather_station)[start_index:end_index]
    weather_df = load_CARRA_data_opt(weather_station)[55200:58199]   
    #weather_df = weather_df.set_index("time")
    #weather_df.index = pd.to_datetime(weather_df.index)
    #time_CARRA = weather_df.index.values
    #print("CARRA, start date: ",weather_df['time'][0], ", end date: ",weather_df['time'][-1] )

    T_CARRA = weather_df.AirTemperature1C.values + 273.15
    z_T_CARRA = weather_df.HeightTemperature1m.values
    RH_CARRA = weather_df.RelativeHumidity1.values
    z_RH_CARRA = weather_df.HeightHumidity1m.values
    WS_CARRA = weather_df.WindSpeed1ms.values
    z_WS_CARRA = weather_df.HeightWindSpeed1m.values
    pres_CARRA = weather_df.AirPressurehPa.values
    SRin_CARRA = weather_df.ShortwaveRadiationDownWm2.values
    SRout_CARRA = weather_df.ShortwaveRadiationUpWm2.values
    LRin_CARRA = weather_df.LongwaveRadiationDownWm2.values
    LRout_CARRA = weather_df.LongwaveRadiationUpWm2.values
    snowfall_CARRA = weather_df.Snowfallmweq.values
    rainfall_CARRA = weather_df.Rainfallmweq.values

    i_AWS = 0
    for i in range(0,10):
       print("compare aws and carra:")
       print("time: ",time_AWS3[i], weather_df.index[i])
       print("temp: ",T_AWS3[i], T_CARRA[i])
       print("z t: ",z_T_AWS3[i], z_T_CARRA[i])
       print("rh: ",RH_AWS3[i], RH_CARRA[i])
       print("rh z: ",z_RH_AWS3[i], z_RH_CARRA[i])
       print("ws: ",WS_AWS3[i], WS_CARRA[i])
       print("z ws: ",z_WS_AWS3[i], z_WS_CARRA[i])
       print("pres: ", pres_AWS3[i], pres_CARRA[i])
       print("sr in: ", SRin_AWS3[i], SRin_CARRA[i])
       print("sr out: ", SRout_AWS3[i], SRout_CARRA[i])
       print("lr in: ", LRin_AWS3[i], LRin_CARRA[i])
       print("lr out: ", LRout_AWS3[i], LRout_CARRA[i])
       print("snow: ", snowfall_AWS3[i], snowfall_CARRA[i])
       print("rain: ", rainfall_AWS3[i], rainfall_CARRA[i])
    

if __name__ == "__main__":
    compare_data()