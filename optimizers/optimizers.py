import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

from cpymad.madx import Madx
from element_parser.data_parser import describe_elements, describe_correctors


class GaussNewton:
    def __init__(self, structure, step, iteration=2):
        self.structure = structure
        self.elements_to_vary, self.initial_parameters = describe_elements(self.structure.structure, "madx\elements\elems.txt")
        self.correctors, _ = describe_correctors(self.structure.structure, "madx\correctors\correctors.txt")
        self.bad_correctors, _ = describe_correctors(self.structure.bad_structure, "madx\correctors\correctors.txt")
        self.bad_elements_to_vary, self.bad_initial_parameters = describe_elements(self.structure.bad_structure, "madx\elements\elems.txt")
        self.elements_number = len(self.elements_to_vary)
        self.shape = [22, self.elements_number]
        # self.shape = [22 * 108, self.elements_number]
        self.step = step
        self.iteration = iteration

    def optimize_lattice(self):
        """
        Calculate necessary parameters changes to make structure correction.

        :param float step: step to vary elements parameters
        :return:
        """
        accumulative_param_additive = np.zeros(self.elements_number)

        bad_response_matrix = self.structure.calculate_response_matrix(self.structure.bad_structure,
                                                                       self.bad_elements_to_vary,
                                                                       accumulative_param_additive,
                                                                       self.bad_correctors)
        model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                         self.elements_to_vary,
                                                                         accumulative_param_additive,
                                                                         self.correctors)
        # initial_vector = np.sum(bad_response_matrix-model_response_matrix, 1)
        initial_vector = np.sum(bad_response_matrix - model_response_matrix, 0)
        # initial_vector = (bad_response_matrix - model_response_matrix).to_numpy().flatten()
        initial_residual = np.sum(initial_vector*initial_vector)

        count = 1
        while count <= self.iteration:
            model_response_matrix_1 = self.structure.calculate_response_matrix(self.structure.structure,
                                                                               self.elements_to_vary,
                                                                               accumulative_param_additive,
                                                                               self.correctors)
            # vector_1 = np.sum(bad_response_matrix-model_response_matrix_1, 1)
            vector_1 = np.sum(bad_response_matrix - model_response_matrix_1, 0)
            # vector_1 = (bad_response_matrix - model_response_matrix_1).to_numpy().flatten()

            J = self.calculate_jacobian(accumulative_param_additive, model_response_matrix_1)
            u, sv, v = self.drop_bad_singular_values(J)
            delta = self.calculate_parameters_delta(J, u, sv, v, vector_1)
            accumulative_param_additive += delta
            count += 1

            fitted_model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                                    self.elements_to_vary,
                                                                                    accumulative_param_additive,
                                                                                    self.correctors)
            # final_vector = np.sum(bad_response_matrix-fitted_model_response_matrix, 1)
            final_vector = np.sum(bad_response_matrix - fitted_model_response_matrix, 0)
            # final_vector = (bad_response_matrix - fitted_model_response_matrix).to_numpy().flatten()
            final_residual = np.sum(final_vector*final_vector)

            print("Initial parameters: ", self.initial_parameters)
            print("Bad parameters: ", self.bad_initial_parameters)
            print("Final parameters: ", list(self.initial_parameters+accumulative_param_additive))
            print("Initial residual: ", initial_residual)
            print("Final residual: ", final_residual)

        return accumulative_param_additive

    def calculate_jacobian(self, accumulative_param_additive, model_response_matrix_1):
        J = np.zeros(self.shape)
        k = 0
        for i in range(self.elements_number):
            now = datetime.now()
            print('Calc Jacob ', str(i) + '/' + str(self.elements_number))
            accumulative_param_variation = np.zeros(self.elements_number)
            accumulative_param_variation[i] = self.step
            accumulative_param_variation += accumulative_param_additive

            model_response_matrix_2 = self.structure.calculate_response_matrix(self.structure.structure,
                                                                               self.elements_to_vary,
                                                                               accumulative_param_variation,
                                                                               self.correctors)
            # vector_2 = np.sum(model_response_matrix_1-model_response_matrix_2, 1)
            vector_2 = np.sum(model_response_matrix_1 - model_response_matrix_2, 0)
            # vector_2 = (model_response_matrix_1 - model_response_matrix_2).to_numpy().flatten()
            J[:,i] = vector_2 / self.step
            print(datetime.now()-now)
            k += 1

        return J

    def drop_bad_singular_values(self, J):
        svd = np.linalg.svd(np.matmul(J.T,J), full_matrices=False)
        u, sv, v = svd
        print("Singulars: ", sv)
        # plt.plot(sv,marker='o')
        # plt.show()
        sv = np.linalg.pinv(np.diag(sv))
        for i in range(len(sv)):
            if sv[i,i] > np.inf:
                sv[i,i] = 0

        return u, sv, v

    def calculate_parameters_delta(self, J, u, sv, v, vector_1):
        J_new = np.matmul(np.matmul(v.T,sv),u.T)
        # delta = -np.matmul(np.linalg.pinv(np.matmul(J.T,J)),J.T).dot(vector_1)
        delta = -np.matmul(J_new,J.T).dot(vector_1)
        # delta = np.where(np.abs(delta)>1e-3, 0.5*delta, delta) # To remove twiss instabillity
        print("delta", delta)

        return delta

    def optimize_orbit(self):
        experimental_responses = self.structure.collect_responses_on_quadrupole_kicks(self.structure.structure,
                                                                                      self.elements_to_vary,
                                                                                      gradient_step=1e-5)

        quadrupole_alignments = {}
        for i, quadrupole in enumerate(self.elements_to_vary):
            response_matrix, response_matrix_at_beginning = self.structure.calculate_response_matrix_on_quadrupole_variation(self.structure.structure,
                                                                                                                             quadrupole,
                                                                                                                             1e-5,
                                                                                                                             i)
            delta = np.linalg.pinv(response_matrix).dot(experimental_responses[quadrupole[0]])
            quadrupole_alignments[quadrupole[0]] = delta
            print('Quadrupole alignment:', str(i) + '/' + str(self.elements_number))
            print('asdd',response_matrix_at_beginning)
            x = response_matrix_at_beginning.dot(delta)
            print("final", x)

        quadrupole_alignments = pd.DataFrame(quadrupole_alignments, index=['dx', 'dy'])

        return quadrupole_alignments


class LevenbergMarquardt(GaussNewton):
    def __init__(self, structure, step, iteration=1, coefficient_lambda=0.0001):
        super().__init__(structure, step, iteration)
        self.coefficient_lambda = coefficient_lambda

    def drop_bad_singular_values(self, J):
        svd = np.linalg.svd(np.matmul(J.T,J)+self.coefficient_lambda*np.diag(np.matmul(J.T,J)), full_matrices=False)
        u, sv, v = svd
        print("Singulars: ", sv)
        sv = np.linalg.pinv(np.diag(sv))
        for i in range(len(sv)):
            if sv[i,i] > np.inf:
                sv[i,i] = 0

        return u, sv, v


class GaussNewtonConstrained(GaussNewton):
    def __init__(self, structure, step, iteration=1):
        super().__init__(structure, step, iteration)
        self.weights = np.ones(self.elements_number)*1e-3

    def optimize_lattice(self):
        """
        Calculate necessary parameters changes to make structure correction.

        :param float step: step to vary elements parameters
        :return:
        """
        accumulative_param_additive = delta = np.zeros(self.elements_number)

        bad_response_matrix = self.structure.calculate_response_matrix(self.structure.bad_structure,
                                                                       self.bad_elements_to_vary,
                                                                       accumulative_param_additive,
                                                                       self.bad_correctors)
        model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                         self.elements_to_vary,
                                                                         accumulative_param_additive,
                                                                         self.correctors)
        # initial_vector = np.sum(bad_response_matrix-model_response_matrix, 1)
        # initial_vector = np.sum(bad_response_matrix - model_response_matrix, 0)
        initial_vector = (bad_response_matrix - model_response_matrix).to_numpy().flatten()
        initial_vector = np.concatenate((initial_vector, self.weights*delta))
        initial_residual = np.sum(initial_vector * initial_vector)

        count = 1
        while count <= self.iteration:
            model_response_matrix_1 = self.structure.calculate_response_matrix(self.structure.structure,
                                                                               self.elements_to_vary,
                                                                               accumulative_param_additive,
                                                                               self.correctors)
            # vector_1 = np.sum(bad_response_matrix-model_response_matrix_1, 1)
            # vector_1 = np.sum(bad_response_matrix - model_response_matrix_1, 0)
            vector_1 = (bad_response_matrix - model_response_matrix_1).to_numpy().flatten()
            vector_1 = np.concatenate((vector_1, self.weights * delta))

            J = self.calculate_jacobian(accumulative_param_additive, model_response_matrix_1)
            J = self.create_modified_jacobian(J, delta)
            u, sv, v = self.drop_bad_singular_values(J)
            delta = self.calculate_parameters_delta(J, u, sv, v, vector_1)
            accumulative_param_additive += delta
            count += 1

            fitted_model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                                    self.elements_to_vary,
                                                                                    accumulative_param_additive,
                                                                                    self.correctors)
            # final_vector = np.sum(bad_response_matrix-fitted_model_response_matrix, 1)
            # final_vector = np.sum(bad_response_matrix - fitted_model_response_matrix, 0)
            final_vector = (bad_response_matrix - fitted_model_response_matrix).to_numpy().flatten()
            final_vector = np.concatenate((final_vector, self.weights * delta))
            final_residual = np.sum(final_vector * final_vector)

            print("Initial parameters: ", self.initial_parameters)
            print("Bad parameters: ", self.bad_initial_parameters)
            print("Final parameters: ", list(self.initial_parameters + accumulative_param_additive))
            print("Initial residual: ", initial_residual)
            print("Final residual: ", final_residual)

        return accumulative_param_additive

    def create_modified_jacobian(self, J, delta_qrad):
        W = np.diag(self.weights * delta_qrad)
        return np.concatenate((J, W), axis=0)
