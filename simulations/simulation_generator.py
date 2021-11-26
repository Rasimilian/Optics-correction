import numpy as np
import json
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

from cpymad.madx import TwissFailed
from madx.madx_tool import Structure
from optimizers.optimizers import GaussNewton, LevenbergMarquardt, GaussNewtonConstrained
from optimizers.py_optimizers import LeastSquaresSolver


def save_plots(structure, optimizer, parameters_delta, path):
    model_twiss_table = structure.twiss_table_4D
    bad_twiss_table = structure.bad_twiss_table_4D
    # model_twiss_table = structure.twiss_table_6D
    # bad_twiss_table = structure.bad_twiss_table_6D
    try:
        corrected_twiss_table = structure.change_structure_for_correction(structure.bad_structure,
                                                                          optimizer.bad_elements_to_vary, -parameters_delta)
        # betx
        plt.plot(model_twiss_table.s, corrected_twiss_table.betx, 'r', label='Corrected')
        plt.plot(model_twiss_table.s, model_twiss_table.betx, 'b', label='Model')
        plt.plot(model_twiss_table.s, bad_twiss_table.betx, linestyle='dashed', color='k', label='Real')
        # plt.plot(model_twiss_table.s, corrected_twiss_table.beta11, 'r', label='Corrected')
        # plt.plot(model_twiss_table.s,model_twiss_table.beta11, 'b', label='Model')
        # plt.plot(model_twiss_table.s,bad_twiss_table.beta11,linestyle='dashed', color='k', label='Real')
        plt.legend()
        plt.xlabel('s, m')
        plt.ylabel('betx, m')
        plt.savefig(path / "betx")
        plt.close()

        # bety
        plt.plot(model_twiss_table.s, corrected_twiss_table.bety, 'r', label='Corrected')
        plt.plot(model_twiss_table.s, model_twiss_table.bety, 'b', label='Model')
        plt.plot(model_twiss_table.s, bad_twiss_table.bety, linestyle='dashed', color='k', label='Real')
        # plt.plot(model_twiss_table.s, corrected_twiss_table.beta22, 'r', label='Corrected')
        # plt.plot(model_twiss_table.s,model_twiss_table.beta22, 'b', label='Model')
        # plt.plot(model_twiss_table.s,bad_twiss_table.beta22,linestyle='dashed', color='k', label='Real')
        plt.legend()
        plt.xlabel('s, m')
        plt.ylabel('bety, m')
        plt.savefig(path / "bety")
        plt.close()

        # beating betx
        plt.plot(model_twiss_table.s, (corrected_twiss_table.betx - model_twiss_table.betx) / model_twiss_table.betx,
                 'r',
                 label='Corrected')
        plt.plot(model_twiss_table.s, (bad_twiss_table.betx - model_twiss_table.betx) / model_twiss_table.betx,
                 linestyle='dashed', color='k', label='Real')
        # plt.plot(model_twiss_table.s,(corrected_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11, 'r', label='Corrected')
        # plt.plot(model_twiss_table.s,(bad_twiss_table.beta11-model_twiss_table.beta11)/model_twiss_table.beta11,linestyle='dashed', color='k', label='Real')

        plt.legend()
        plt.xlabel('s, m')
        plt.ylabel('x beta-beating, m')
        plt.savefig(path / "beating_betx")
        plt.close()

        # beating bety
        plt.plot(model_twiss_table.s, (corrected_twiss_table.bety - model_twiss_table.bety) / model_twiss_table.bety,
                 'r',
                 label='Corrected')
        plt.plot(model_twiss_table.s, (bad_twiss_table.bety - model_twiss_table.bety) / model_twiss_table.bety,
                 linestyle='dashed', color='k', label='Real')
        # plt.plot(model_twiss_table.s,(corrected_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22, 'r', label='Corrected')
        # plt.plot(model_twiss_table.s,(bad_twiss_table.beta22-model_twiss_table.beta22)/model_twiss_table.beta22,linestyle='dashed', color='k', label='Real')

        plt.legend()
        plt.xlabel('s, m')
        plt.ylabel('y beta-beating, m')
        plt.savefig(path / "beating_bety")
        plt.close()
    except TwissFailed:
        print("Cant make plots")






def make_simulation(output_dir, algorithm):
    root_lat = "madx/structures"
    root_elem = "madx/elements"
    ideal_lattice = "VEPP4M_full1.txt"
    erroneous_lattices = ["VEPP4M_full1_3_quads_errors.txt",
                          "VEPP4M_full1_combined_magnets_errors.txt",
                          "VEPP4M_full1_all_grads_errors.txt"]
    erroneous_lattices = ["VEPP4M_full1_all_grads_errors.txt"]
    element_sets = ["combined_magnets.txt",
                    "elems.txt",
                    "quads.txt"]
    element_sets = ["elems.txt"]
    # grad_steps = [1e-3, 1e-4]
    grad_steps = [1e-3]
    # corrector_steps = [1e-3, 1e-4, 1e-5, 1e-6]
    corrector_steps = [1e-6]
    matrix_forms = ["Coupled", "Uncoupled"]
    matrix_forms = ["Coupled"]

    if algorithm == "GaussNewton":
        optimizer = GaussNewton
        additional_parameters = [None]
    elif algorithm == "LevenbergMarquardt":
        optimizer = LevenbergMarquardt
        # additional_parameters = [1e-1, 1e-2, 1e-3, 1e-4]  # lambda param
        additional_parameters = [1e-3]
    else:
        optimizer = GaussNewtonConstrained
        additional_parameters = [1e-1, 1e-2, 1e-3, 1e-4]  # weigths

    log = []
    sim_num = 0

    for erroneous_lattice in erroneous_lattices:
        for element_set in element_sets:
            for grad_step in grad_steps:
                for corrector_step in corrector_steps:
                    for additional_parameter in additional_parameters:
                        for matrix_form in matrix_forms:
                            sim_num += 1
                            print("sim_num:", sim_num)
                            start_time = datetime.now()
                            try:
                                structure = Structure(structure_file=root_lat + "/" + ideal_lattice,
                                                      bad_structure_file=root_lat + "/" + erroneous_lattice)
                                solver = optimizer(structure, grad_step, corrector_step, 3,
                                                   root_elem + "/" + element_set, additional_parameter)
                                parameters_delta, history = solver.optimize_lattice(matrix_form)
                                end_time = datetime.now()

                                res_to_dict = {"Lattice": erroneous_lattice,
                                               "Elements": element_set,
                                               "Optimizer": algorithm,
                                               "Gradient step": grad_step,
                                               "Corrector step": corrector_step,
                                               "Weigths_or_lambdas": additional_parameter,
                                               "Matrix form": matrix_form,
                                               "Execution time": str(end_time - start_time),
                                               "Results": history
                                               }
                                log.append(res_to_dict)
                                print(end_time)

                                output_dir_final = output_dir / f"sim_{sim_num}_{Path(erroneous_lattice).stem}_{Path(element_set).stem}_{algorithm}"
                                output_dir_final.mkdir(parents=True, exist_ok=True)
                                save_plots(structure, solver, parameters_delta, output_dir_final)
                                with open(output_dir_final / "results.json", 'w') as f:
                                    json.dump(res_to_dict, f, indent=4)
                            except TwissFailed:
                                print("Cant make sim")



if __name__ == "__main__":
    output_dir = Path("simulations/results") / str(datetime.now().strftime("%d_%m_%y-%H-%M-%S"))
    # make_simulation(output_dir, "GaussNewton")
    # algorithms = ["GaussNewton", "LevenbergMarquardt", "GaussNewtonConstrained"]
    algorithms = ["LevenbergMarquardt"]
    for algorithm in algorithms:
        make_simulation(output_dir, algorithm)
