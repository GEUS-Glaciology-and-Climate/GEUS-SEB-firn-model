import __init__
import time
import main_SEB_firn as main
import lib_seb_smb_model as seb
import numpy as np

def run_measure():
    time_array = []

    for i in range(30):
        try:
            cpu_start_time = time.process_time()
            main.run_SEB_firn()
            cpu_end_time = time.process_time()
            cpu_time = (cpu_end_time - cpu_start_time)
            time_array.append(cpu_time)
            print(time_array)
        except:
            print ("error "+str(IOError))
            pass

    for i in range(len(time_array)):
        print(time_array[i])

    print(time_array)
    print("Mean:")
    print(str(np.mean(time_array)))

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
    c = main.setConstants()

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

    print("Optimized code: ")
   # print(time_opt)
    print("Mean:")
    print(str(np.mean(time_opt)))

    print("Not optimized code: ")
   # print(time_not_opt)
    print("Mean:")
    print(str(np.mean(time_not_opt)))


measure_SensLatFluxes()

#run_measure()



