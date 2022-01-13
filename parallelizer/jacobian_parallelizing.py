import ctypes
import numpy as np
import multiprocessing as mp
from datetime import datetime
from cpymad.madx import Madx


def add_worker_process(calc_jacob, model_response_matrix, accumulative_param_additive, queue, shape, arrD, structure, step, elements_to_vary, correctors):
    J = np.reshape(np.frombuffer(arrD), shape)

    while True:
        # job is an array of params portion (indices of params)
        job = queue.get()
        if not isinstance(job, np.ndarray):
            break

        J[:, job] = calc_jacob(accumulative_param_additive, model_response_matrix, job, structure, step, elements_to_vary, correctors)

        queue.task_done()
    queue.task_done()


def parallelize(calc_jacob, model_response_matrix, accumulative_param_additive, shape, structure, step, elements_to_vary, correctors):
    arrD = mp.RawArray(ctypes.c_double, shape[0] * shape[1])

    # setup jobs
    # nCPU = mp.cpu_count()
    nCPU = 4
    nJobs = nCPU

    param_portions = np.array_split(np.arange(len(accumulative_param_additive)), nJobs)

    queue = mp.JoinableQueue()
    for param_portion in param_portions:
        queue.put(param_portion)
    for i in range(nCPU):
        queue.put(None)

    # run workers
    workers = []
    for i in range(nCPU):
        worker = mp.Process(target=add_worker_process,
                            args=(calc_jacob, model_response_matrix, accumulative_param_additive, queue, shape, arrD, structure, step, elements_to_vary, correctors))
        workers.append(worker)
        worker.start()

    queue.join()

    # make array from shared memory
    J = np.reshape(np.frombuffer(arrD), shape)

    return J
