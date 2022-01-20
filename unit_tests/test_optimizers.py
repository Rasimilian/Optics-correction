import pytest
from madx.madx_tool import Structure
from optimizers.optimizers import LevenbergMarquardt
from numpy.testing import assert_array_equal, assert_array_almost_equal


@pytest.mark.parametrize(
    'structure_file, bad_structure_file, correctors, elements_to_vary, iteration, corrector_step, grad_step, coefficient_lambda, expected_results',
    [
        ('madx/structures/VEPP4M_full1.txt', 'madx/structures/VEPP4M_full1_3_quads_errors.txt',
         'madx/correctors/correctors.txt',
         'madx/elements/quads.txt', 3, 1e-6, 1e-3, 1e-3,
         ([
              0.23857, -0.27657, 0.36024, -0.31948, 0.93011966, -0.8932828, 0.99208665, -0.8760174, -0.27743, 0.1809,
              -0.24116, 0.34942, -0.69308, -0.69308, 0.34942, -0.24116, 0.1809, -0.27743, -0.8760174, 0.99208665,
              -0.8932828, 0.93011966, -0.31948, 0.36024, -0.27657, 0.23857
          ],
          [
              0.24157, -0.27657, 0.36024, -0.31948, 0.94011966, -0.8832828, 0.99208665, -0.8760174, -0.27743, 0.1809,
              -0.24116, 0.34942, -0.69308, -0.69308, 0.34942, -0.24116, 0.1809, -0.27743, -0.8760174, 0.99208665,
              -0.8932828, 0.93011966, -0.31948, 0.36024, -0.27657, 0.23857
          ],
          [
              0.003, 0.0, 0.0, 0.0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0
          ],
          18.478947471750622, 5.147806556337487e-09)
         )
    ])
def test_levenberg_marquardt(structure_file,
                             bad_structure_file,
                             correctors,
                             elements_to_vary,
                             iteration,
                             corrector_step,
                             grad_step,
                             coefficient_lambda,
                             expected_results):
    structure = Structure(structure_file, bad_structure_file)
    optimizer = LevenbergMarquardt(structure, correctors, elements_to_vary, corrector_step, grad_step, iteration,
                                   coefficient_lambda)
    parameters_delta, initial_residual, final_residual = optimizer.optimize_lattice()

    assert_array_equal(optimizer.initial_parameters, expected_results[0])
    assert_array_equal(optimizer.bad_initial_parameters, expected_results[1])
    assert_array_almost_equal(parameters_delta, expected_results[2], 4)
    assert initial_residual == expected_results[3]
    assert final_residual == expected_results[4]

    # model_twiss_table = structure.twiss_table_4D
    # bad_twiss_table = structure.bad_twiss_table_4D
    # corrected_twiss_table = structure.change_structure_for_correction(structure.bad_structure,
    #                                                                   optimizer.bad_elements_to_vary, -parameters_delta)
