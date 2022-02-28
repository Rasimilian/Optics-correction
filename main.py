import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import cProfile

from madx.madx_tool import Structure
from optimizers.optimizers import GaussNewton, LevenbergMarquardt, GaussNewtonConstrained
from optimizers.py_optimizers import LeastSquaresSolver
from multiprocessing import Pool

if __name__ == "__main__":
    structure = Structure()

    now = datetime.now()
    # optimizer = GaussNewton(structure, step=1e-3)
    optimizer = GaussNewton(structure, "madx/correctors/correctors.txt", "madx/elements/quads.txt", 1e-6,
                                   grad_step=1e-3)
    # optimizer = LevenbergMarquardt(structure, "madx/correctors/correctors.txt", "madx/elements/quads.txt", 1e-6, grad_step=1e-3)
    # optimizer = GaussNewtonConstrained(structure, step=3e-3)

    parameters_delta, _, _, alignments_delta = optimizer.optimize_lattice()
    # parameters_delta = optimizer.optimize_orbit()
    print(parameters_delta)
    print(alignments_delta)
    print(datetime.now()-now)

    model_twiss_table = structure.twiss_table_4D
    bad_twiss_table = structure.bad_twiss_table_4D
    # model_twiss_table = structure.twiss_table_6D
    # bad_twiss_table = structure.bad_twiss_table_6D
    corrected_twiss_table = structure.change_structure_for_correction(structure.bad_structure, optimizer.bad_elements_to_vary, -parameters_delta, alignments_delta, optimizer.names)

    plt.plot(model_twiss_table.s,corrected_twiss_table.betx, 'r', label='Corrected')
    plt.plot(model_twiss_table.s,model_twiss_table.betx, 'b', label='Model')
    plt.plot(model_twiss_table.s,bad_twiss_table.betx,linestyle='dashed', color='k', label='Real')
    # plt.plot(model_twiss_table.s,model_twiss_table.beta11, 'b', label='Model')
    # plt.plot(model_twiss_table.s,bad_twiss_table.beta11,linestyle='dashed', color='k', label='Real')
    plt.legend()
    plt.xlabel('s, m')
    plt.ylabel(r'$\beta_x$, m')
    plt.show()

    plt.plot(model_twiss_table.s,corrected_twiss_table.bety, 'r', label='Corrected')
    plt.plot(model_twiss_table.s,model_twiss_table.bety, 'b', label='Model')
    plt.plot(model_twiss_table.s,bad_twiss_table.bety,linestyle='dashed', color='k', label='Real')
    # plt.plot(model_twiss_table.s,model_twiss_table.beta22, 'b', label='Model')
    # plt.plot(model_twiss_table.s,bad_twiss_table.beta22,linestyle='dashed', color='k', label='Real')
    plt.legend()
    plt.xlabel('s, m')
    plt.ylabel(r'$\beta_y$, m')
    plt.show()

    plt.plot(model_twiss_table.s,(corrected_twiss_table.betx-model_twiss_table.betx)/model_twiss_table.betx, 'r', label='Corrected')
    plt.plot(model_twiss_table.s,(bad_twiss_table.betx-model_twiss_table.betx)/model_twiss_table.betx,linestyle='dashed', color='k', label='Real')
    # plt.plot(model_twiss_table.s,(corrected_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11, 'r', label='Corrected')
    # plt.plot(model_twiss_table.s,(bad_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11,linestyle='dashed', color='k', label='Real')
    plt.legend()
    plt.xlabel('s, m')
    plt.ylabel(r'$\Delta\beta_x/\beta_x$, m')
    plt.show()

    plt.plot(model_twiss_table.s,(corrected_twiss_table.bety-model_twiss_table.bety)/model_twiss_table.bety, 'r', label='Corrected')
    plt.plot(model_twiss_table.s,(bad_twiss_table.bety-model_twiss_table.bety)/model_twiss_table.bety,linestyle='dashed', color='k', label='Real')
    # plt.plot(model_twiss_table.s,(corrected_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22, 'r', label='Corrected')
    # plt.plot(model_twiss_table.s,(bad_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22,linestyle='dashed', color='k', label='Real')
    plt.legend()
    plt.xlabel('s, m')
    plt.ylabel(r'$\Delta\beta_y/\beta_y$, m')
    plt.show()
