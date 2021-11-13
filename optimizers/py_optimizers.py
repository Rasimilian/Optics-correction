import numpy as np
from scipy.optimize import least_squares, minimize, minimize_scalar

from cpymad.madx import Madx
from element_parser.data_parser import describe_elements, describe_correctors


class LeastSquaresSolver:
    def __init__(self, structure):
        self.structure = structure
        self.elements_to_vary, self.initial_parameters = describe_elements(self.structure.structure, "madx\elements\elems.txt")
        self.correctors, _ = describe_correctors(self.structure.structure, "madx\correctors\correctors.txt")
        self.bad_correctors, _ = describe_correctors(self.structure.bad_structure, "madx\correctors\correctors.txt")
        self.bad_elements_to_vary, self.bad_initial_parameters = describe_elements(self.structure.bad_structure, "madx\elements\elems.txt")
        self.elements_number = len(self.elements_to_vary)
        # self.shape = [22, self.elements_number]
        self.shape = [22 * 108, self.elements_number]

    def get_target_function(self, params, bad_response_matrix):
        model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                         self.elements_to_vary,
                                                                         params,
                                                                         self.correctors)
        # vector = np.sum(bad_response_matrix - model_response_matrix, 0)
        vector = (bad_response_matrix - model_response_matrix).to_numpy().flatten()
        residual = np.sum(vector * vector)
        return residual

    def optimize(self):
        bad_response_matrix = self.structure.calculate_response_matrix(self.structure.bad_structure,
                                                                       self.bad_elements_to_vary,
                                                                       np.zeros(self.elements_number),
                                                                       self.bad_correctors)
        initial_guess = [x - y for x, y in zip(self.bad_initial_parameters, self.initial_parameters)]
        print(initial_guess)
        initial_guess = np.zeros(self.elements_number)
        res = least_squares(self.get_target_function,
                            initial_guess,
                            args=(bad_response_matrix,),
                            verbose=2,
                            ftol=1e-8,
                            x_scale='jac',
                            method='dogbox')
        print(res)

    def optimize1(self):
        bad_response_matrix = self.structure.calculate_response_matrix(self.structure.bad_structure,
                                                                       self.bad_elements_to_vary,
                                                                       np.zeros(self.elements_number),
                                                                       self.bad_correctors)
        initial_guess = [x - y for x, y in zip(self.bad_initial_parameters, self.initial_parameters)]
        print(initial_guess)
        initial_guess = np.zeros(self.elements_number)
        res = minimize(self.get_target_function,
                       initial_guess,
                       args=(bad_response_matrix,),
                       method='Nelder-Mead')
        print(res)

    def optimize2(self):
        bad_response_matrix = self.structure.calculate_response_matrix(self.structure.bad_structure,
                                                                       self.bad_elements_to_vary,
                                                                       np.zeros(self.elements_number),
                                                                       self.bad_correctors)
        initial_guess = [x - y for x, y in zip(self.bad_initial_parameters, self.initial_parameters)]
        print(initial_guess)
        initial_guess = np.zeros(self.elements_number)
        res = minimize_scalar(self.get_target_function,
                              initial_guess,
                              args=(bad_response_matrix,))
        print(res)
