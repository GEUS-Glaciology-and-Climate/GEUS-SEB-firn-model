import __init__
import time
#from main_SEB_firn import run_SEB_firn
import main_SEB_firn as main
import numpy as np

# arr = np.array([54.640625, 54.796875, 54.84375, 55.296875, 55.078125, 55.234375, 54.953125, 54.90625, 55.3125, 54.625])
# print(arr)
# print(np.mean(arr))

def run_measure():
    time_array = []

    for i in range(5):
        cpu_start_time = time.process_time()
        main.run_SEB_firn()
        cpu_end_time = time.process_time()
        cpu_time = (cpu_end_time - cpu_start_time)
        time_array.append(cpu_time)
        print(time_array)

    print(time_array)
    print("Mean: "+ str(np.mean(time_array)))

run_measure()


