import __init__
import lib_io as io
import lib_seb_smb_model as seb
import numpy as np
import pandas as pd
import time

from main_SEB_firn import run_SEB_firn, set_constants
from main_firn import run_main_firn_old, run_main_firn_parallel

def run_measure_SEB_firn():
    time_array = []
    print("Measuring run_SEB_firn (main_SEB_firn). \nHave you controlled num of layers? To plot or not to plot?")

    for i in range(10):
        # Run test, keep going if error
        try:
            cpu_start_time = time.process_time()
            run_SEB_firn()
            cpu_end_time = time.process_time()
            cpu_time = (cpu_end_time - cpu_start_time)
            time_array.append(cpu_time)
            print("time in measure: ", time_array)
        except Exception as e:
            print (e)
            pass

        # Run test, breaks if error
        # cpu_start_time = time.process_time()
        # run_SEB_firn()
        # cpu_end_time = time.process_time()
        # cpu_time = (cpu_end_time - cpu_start_time)
        # time_array.append(cpu_time)
        # print(time_array)

    # for i in range(len(time_array)):
    #     print(time_array[i])

    print(time_array)
    print("Mean:")
    print(str(np.mean(time_array)))


def measure_firn_parallel():
    print("--- Measure firn parallel running ---")
    site_list = ["KAN_U", "KAN_M"]
    time_list = []
    numb_loops = 10
    
    for i in range(numb_loops):
        #Parallel running
        parallel_start_time = time.process_time()
        run_main_firn_parallel(site_list)
        parallel_end_time = time.process_time()
        parallel_time = (parallel_end_time - parallel_start_time)
        time_list.append(parallel_time)

    return time_list


def measure_firn_old():
    print("Are you plotting anything?")
    print("--- Measure firn old running ---")
    site_list = ["KAN_U", "KAN_M", "KAN_U", "KAN_M"]

    time_list = []
    tot_time = 0
    numb_loops = 10
    for i in range(numb_loops):
        time_old_start = time.process_time()
        run_main_firn_old(site_list)
        time_old_end = time.process_time()
        old_time = (time_old_end - time_old_start)
        time_list.append(old_time)
        tot_time = old_time + tot_time

    mean = tot_time / numb_loops 
    std = np.std(np.array(time_list))

    print("Old code, mean and std:")
    print(mean)
    print(std)
    return(time_list)

def run_all_measures_firn():
    # old_loop_time = measure_firn_old()
    # print("Old loop time:")
    # print(old_loop_time)

    parallel_time = measure_firn_parallel()
    print("Parallel time: ")
    print(parallel_time)


def measure_SensLatFluxes():
    '''Measure time diff between optimized and not optimized SensLatFluxes,
      looping 500 times, with synthetic values'''

    time_not_opt = []
    time_opt = []
    WS= 4.88
    nu= 1.5471445341359545e-05
    q = 0.002983881449656623
    snowthick = 1.0
    Tsurf = 50.00032543036966
    theta = 269.23442297512435
    theta_v = 269.72264072580384
    pres = 857.105524
    rho_atm = 1.109137923478789
    z_WS = 2.9995
    z_T = 2.4995
    z_RH = 2.4995
    z_0 = 0.0013
    c = set_constants("KAN_U")

    for i in range(1000):
        # Measure with optimization
        cpu_start_opt = time.process_time()
        for i in range(300):
            seb.SensLatFluxes_bulk_opt(
                    WS,
                    nu,
                    q,
                    snowthick,
                    Tsurf,
                    theta,
                    theta_v,
                    pres,
                    rho_atm,
                    z_WS,
                    z_T,
                    z_RH,
                    z_0,
                    c
                )
        cpu_end_opt = time.process_time()
        cpu_time_opt = (cpu_end_opt - cpu_start_opt)
        time_opt.append(cpu_time_opt)

        # Measure without optimization
        cpu_start_notopt = time.process_time()
        for i in range(300):
            seb.SensLatFluxes_bulk_old(
                    WS,
                    nu,
                    q,
                    snowthick,
                    Tsurf,
                    theta,
                    theta_v,
                    pres,
                    rho_atm,
                    z_WS,
                    z_T,
                    z_RH,
                    z_0,
                    c
                )
        cpu_end_notopt = time.process_time()
        cpu_time_notopt = (cpu_end_notopt - cpu_start_notopt)
        time_not_opt.append(cpu_time_notopt)

    print("Optimized code, mean:")
    print(str(np.mean(time_opt)))

    print("Not optimized code, mean:")
    print(str(np.mean(time_not_opt)))


if __name__ == '__main__':
    #measure_SensLatFluxes()
    #run_measure_SEB_firn()
    #run_measure_firn()
    #compare_full_simulation()
    #check_diff()
    run_all_measures_firn()
    

