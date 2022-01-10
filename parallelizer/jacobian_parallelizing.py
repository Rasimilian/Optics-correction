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

    params_in_job = int(len(accumulative_param_additive) / nJobs)
    # param_portions = np.array_split(np.arange(len(accumulative_param_additive)), params_in_job)
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

# def add_worker_process(calc_jacob, accumulative_param_additive, queue, shape, arrD):
#     J = np.reshape(np.frombuffer(arrD), shape)
#
#     while True:
#         # job is an array of params portion (indices of params)
#         job = queue.get()
#         if not isinstance(job, np.ndarray):
#             break
#
#         for j in job:
#             J[:, j] = calc_jacob(accumulative_param_additive)
#
#         queue.task_done()
#     queue.task_done()
#
#
# def parallelize(calc_jacob, accumulative_param_additive, shape):
#     arrD = mp.RawArray(ctypes.c_double, shape[0] * shape[1])
#
#     # setup jobs
#     # nCPU = mp.cpu_count()
#     nCPU = 2
#     nJobs = nCPU * 5
#
#     params_in_job = int(len(accumulative_param_additive) / nJobs)
#     param_portions = np.array_split(accumulative_param_additive, params_in_job)
#
#     queue = mp.JoinableQueue()
#     for param_portion in param_portions:
#         queue.put(param_portion)
#     for i in range(nCPU):
#         queue.put(None)
#
#     # run workers
#     workers = []
#     for i in range(nCPU):
#         worker = mp.Process(target=add_worker_process,
#                             args=(calc_jacob, accumulative_param_additive, queue, shape, arrD))
#         workers.append(worker)
#         worker.start()
#
#     queue.join()
#
#     # make array from shared memory
#     J = np.reshape(np.frombuffer(arrD), shape)
#     return J

# def func(arr):
#     sum_ = 0
#     for i in arr:
#         sum_ += i
#         print(i)
#     print(sum_)
#     return sum_
#
# def calc_jacobian(accumulative_param_additive):
#     madx = Madx(stdout=False)
#     madx.option(echo=False, warn=False, info=False, twiss_print=False)
#     madx.call(file="madx\structures\VEPP4M_full1.txt")
#     madx.input('beam,particle=electron,energy=1.8;')
#     madx.input('use,sequence=RING;')
#
#     # TODO add errors
#     # self.errors_table = self.add_errors_to_structure(madx, 'quadrupole', value=0.000001, initial_imperfections=initial_imperfections)
#
#     # madx.input('select,flag=twiss,class=monitor;')
#
#     madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
#     madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_BPMs;')
#     twiss_table = madx.table.twiss_in_BPMs.betx
#     madx.quit()
#     return twiss_table



if __name__ == "__main__":
    # arr = np.arange(1, 10000001)
    # parts = np.array_split(arr, 2)
    # time1 = datetime.now()
    # p1 = mp.Process(target=func, args=(parts[0],))
    # p2 = mp.Process(target=func, args=(parts[1],))
    # p1.start()
    # p2.start()
    #
    # p1.join()
    # p2.join()
    #
    # time2 = datetime.now() - time1
    #
    # time3 = datetime.now()
    # # ss = func(arr)
    # time4 = datetime.now() - time3
    #
    # print(time2)
    # print(time4)
    accumulative_param_additive = np.arange(1,100)
    J = parallelize(calc_jacobian, accumulative_param_additive, shape=(651, 100))
    print(J)

