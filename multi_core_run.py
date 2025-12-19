# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

Designed to be run on geussmb01
"""

import os
import multiprocessing


def run_on_core(station, core):
    """Set CPU affinity and run run_SEB_firn for the given station."""
    print(f"Core {core} processing station: {station}")
    surface_input_path = f"/data/CARRA/extracted/list_pixels_minimal_grid/{station}.nc"
    output_path = '/data/CARRA-SMB/list_pixels_minimal_grid/'
    spin_up_path = '/data/CARRA-SMB/spin up 3H/'


    os.system(
        f"taskset -c {core} python3 -c 'import main_SEB_firn; "
        f"main_SEB_firn.run_SEB_firn(\"{station}\, "
        f"\"{surface_input_path}\", \"{output_path}\", \"{spin_up_path}\")'"
    )
    print(f"Core {core} finished station: {station}")

def worker(task_queue):
    """Worker process that executes tasks sequentially on its assigned core."""
    core = task_queue.get()  # Get the assigned core
    while True:
        station = task_queue.get()  # Get the next station to process
        if station is None:
            break  # Stop worker when None is received
        run_on_core(station, core)

def standard_run_parallel(station_list):
    # max_core_usage =  multiprocessing.cpu_count()-1 # all cores except one
    max_core_usage = 21  # Limit to 6 cores
    num_cores = min(len(station_list), max_core_usage)
    task_queues = [multiprocessing.Queue() for _ in range(num_cores)]
    processes = []

    # Start workers with dedicated cores
    for core, task_queue in enumerate(task_queues):
        task_queue.put(core)  # Assign core number to worker
        p = multiprocessing.Process(target=worker, args=(task_queue,))
        p.start()
        processes.append((p, task_queue))

    # Distribute tasks in a round-robin fashion
    for i, station in enumerate(station_list):
        task_queues[i % num_cores].put(station)

    # Send termination signal (None) to workers
    for _, task_queue in processes:
        task_queue.put(None)

    # Wait for all workers to finish
    for p, _ in processes:
        p.join()

if __name__ == "__main__":
    station_list = [s.replace('.nc', '') for s in os.listdir("/data/CARRA/extracted/list_pixels_minimal_grid/")]
    # station_list = [s.replace('.nc', '') for s in os.listdir("./input/weather data/CARRA_at_AWS")]
    # station_list.sort()


    standard_run_parallel(station_list)
