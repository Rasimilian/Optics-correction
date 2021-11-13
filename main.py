import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import cProfile

from madx.madx_tool import Structure
from optimizers.optimizers import GaussNewton, LevenbergMarquardt, GaussNewtonConstrained
from optimizers.py_optimizers import LeastSquaresSolver
from multiprocessing import Pool


structure = Structure()

now = datetime.now()

optimizer = LeastSquaresSolver(structure)
optimizer.optimize1()
print(datetime.now()-now)
#
breakpoint()
# response_matrix = structure.calculate_response_matrix(structure.structure, structure.structure_in_lines, 1e-4, 0)
# print(datetime.now()-now)
# plt.plot(structure.twiss_table_4D.s,structure.twiss_table_4D.x,'r')
# plt.plot(structure.twiss_table_6D.s,structure.twiss_table_6D.x,'b')
# plt.show()

now = datetime.now()
optimizer = GaussNewton(structure, step=1e-3)
# optimizer = LevenbergMarquardt(structure, step=1e-3)
# optimizer = GaussNewtonConstrained(structure, step=1e-3)
parameters_delta = optimizer.optimize_lattice()
# parameters_delta = optimizer.optimize_orbit()
print(parameters_delta)
print(datetime.now()-now)
# breakpoint()
#
#
model_twiss_table = structure.twiss_table_4D
bad_twiss_table = structure.bad_twiss_table_4D
# model_twiss_table = structure.twiss_table_6D
# bad_twiss_table = structure.bad_twiss_table_6D
# corrected_twiss_table = structure.change_structure_for_correction(structure.structure, optimizer.elements_to_vary, parameters_delta)
corrected_twiss_table = structure.change_structure_for_correction(structure.bad_structure, optimizer.bad_elements_to_vary, -parameters_delta)

plt.plot(model_twiss_table.s,corrected_twiss_table.betx, 'r', label='Corrected')
plt.plot(model_twiss_table.s,model_twiss_table.betx, 'b', label='Model')
plt.plot(model_twiss_table.s,bad_twiss_table.betx,linestyle='dashed', color='k', label='Real')
# plt.plot(model_twiss_table.s,model_twiss_table.beta11, 'b', label='Model')
# plt.plot(model_twiss_table.s,bad_twiss_table.beta11,linestyle='dashed', color='k', label='Real')
plt.legend()
plt.xlabel('s, m')
plt.ylabel('betx, m')
plt.show()

plt.plot(model_twiss_table.s,corrected_twiss_table.bety, 'r', label='Corrected')
plt.plot(model_twiss_table.s,model_twiss_table.bety, 'b', label='Model')
plt.plot(model_twiss_table.s,bad_twiss_table.bety,linestyle='dashed', color='k', label='Real')
# plt.plot(model_twiss_table.s,model_twiss_table.beta22, 'b', label='Model')
# plt.plot(model_twiss_table.s,bad_twiss_table.beta22,linestyle='dashed', color='k', label='Real')
plt.legend()
plt.xlabel('s, m')
plt.ylabel('bety, m')
plt.show()

# plt.plot(model_twiss_table.s,model_twiss_table.betx, 'b', label='Model')
plt.plot(model_twiss_table.s,(corrected_twiss_table.betx-model_twiss_table.betx)/model_twiss_table.betx, 'r', label='Corrected')
plt.plot(model_twiss_table.s,(bad_twiss_table.betx-model_twiss_table.betx)/model_twiss_table.betx,linestyle='dashed', color='k', label='Real')
# plt.plot(model_twiss_table.s,(corrected_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11, 'r', label='Corrected')
# plt.plot(model_twiss_table.s,(bad_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11,linestyle='dashed', color='k', label='Real')
plt.legend()
plt.xlabel('s, m')
plt.ylabel('x beta-beating, m')
plt.show()

# plt.plot(model_twiss_table.s,model_twiss_table.betx, 'b', label='Model')
plt.plot(model_twiss_table.s,(corrected_twiss_table.bety-model_twiss_table.bety)/model_twiss_table.bety, 'r', label='Corrected')
plt.plot(model_twiss_table.s,(bad_twiss_table.bety-model_twiss_table.bety)/model_twiss_table.bety,linestyle='dashed', color='k', label='Real')
# plt.plot(model_twiss_table.s,(corrected_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22, 'r', label='Corrected')
# plt.plot(model_twiss_table.s,(bad_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22,linestyle='dashed', color='k', label='Real')
plt.legend()
plt.xlabel('s, m')
plt.ylabel('y beta-beating, m')
plt.show()

print('errors in STL1 0.24157, SIL1 0.94011966, SIL2 -0.8832828')

# def calc(x):
#     # print(x**2)
#     return x**2
#
# if __name__ =='__main__':
#     structure = Structure()
#     structures = [structure.structure,structure.structure,structure.structure]
#     # asd = [1,2,3,4,5]
#     pool = Pool(processes=3)
#     res = pool.map(structure.calculate_structure_4D,structures)
#     # print(pool.map(calc, asd))

