import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from typing import List, Optional, Union, Tuple

from cpymad.madx import Madx, TwissFailed
from madx.madx_tool import Structure
from element_parser.data_parser import describe_elements, describe_correctors, match_elements_indices_of_two_structures


class GaussNewton:
    def __init__(self, structure: Structure, step: float, iteration: int = 2):
        self.structure = structure
        self.elements_to_vary, self.initial_parameters = describe_elements(self.structure.structure,
                                                                           "madx\elements\elems.txt")
        self.correctors, _ = describe_correctors(self.structure.structure, "madx\correctors\correctors.txt")
        self.bad_correctors, _ = describe_correctors(self.structure.bad_structure, "madx\correctors\correctors.txt")
        self.bad_elements_to_vary, self.bad_initial_parameters = describe_elements(self.structure.bad_structure,
                                                                                   "madx\elements\elems.txt")
        self.bad_elements_to_vary, self.bad_initial_parameters = match_elements_indices_of_two_structures(
            self.elements_to_vary, self.initial_parameters, self.bad_elements_to_vary, self.bad_initial_parameters)
        self.elements_number = len(self.elements_to_vary)
        # self.shape = [22, self.elements_number]
        self.shape = [22 * 108, self.elements_number]
        self.step = step
        self.iteration = iteration

    def _get_residual(self,
                      bad_response_matrix: pd.DataFrame,
                      model_response_matrix: pd.DataFrame,
                      weights: Optional[Union[List, np.ndarray]] = None,
                      vector_type: str = "Flatenned") -> Tuple[np.ndarray, float]:
        if vector_type == "axis-0 summed":
            vector = np.sum(bad_response_matrix - model_response_matrix, 0)
        elif vector_type == "axis-1 summed":
            vector = np.sum(bad_response_matrix - model_response_matrix, 1)
        else:
            vector = (bad_response_matrix - model_response_matrix).to_numpy().flatten()

        if isinstance(weights, (list, np.ndarray)):
            vector = np.concatenate((vector, weights))
        residual = np.sum(vector * vector)

        return vector, residual

    def optimize_lattice(self) -> np.ndarray:
        """
        Calculate necessary parameters changes to make structure correction.

        :return: parameters deltas which fit erroneous structure
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
        initial_vector, initial_residual = self._get_residual(bad_response_matrix, model_response_matrix)

        final_residual = 1000
        count = 1
        # while count <= self.iteration:
        while final_residual > 1:
            model_response_matrix_1 = self.structure.calculate_response_matrix(self.structure.structure,
                                                                               self.elements_to_vary,
                                                                               accumulative_param_additive,
                                                                               self.correctors)
            vector_1, _ = self._get_residual(bad_response_matrix, model_response_matrix_1)

            J = self.calculate_jacobian(accumulative_param_additive, model_response_matrix_1)
            u, sv, v = self.drop_bad_singular_values(J)
            delta = self.calculate_parameters_delta(J, u, sv, v, vector_1)
            accumulative_param_additive += delta
            count += 1

            is_fitted = False
            k = 1
            while not is_fitted:
                try:
                    fitted_model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                                            self.elements_to_vary,
                                                                                            accumulative_param_additive,
                                                                                            self.correctors)
                    is_fitted = True
                except TwissFailed:
                    accumulative_param_additive = np.where(np.abs(accumulative_param_additive) < k, accumulative_param_additive, 0)
                    k /= 2
                    print("deltas reduced")
                    print(k)
                    print(accumulative_param_additive)
            final_vector, final_residual = self._get_residual(bad_response_matrix, fitted_model_response_matrix)

            print("Iteration: ", count)
            print("Initial parameters: ", self.initial_parameters)
            print("Bad parameters: ", self.bad_initial_parameters)
            print("Final parameters: ", list(self.initial_parameters + accumulative_param_additive))
            print("Initial residual: ", initial_residual)
            print("Final residual: ", final_residual)

        return accumulative_param_additive

    def calculate_jacobian(self, accumulative_param_additive: np.ndarray,
                           model_response_matrix_1: pd.DataFrame) -> np.ndarray:
        """
        Function to calculate Jacobian.

        :param accumulative_param_additive: parameters deltas from last iteration
        :param model_response_matrix_1: interim response matrix
        :return: Jacobian
        """
        J = np.zeros(self.shape)
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
            vector_2, _ = self._get_residual(model_response_matrix_1, model_response_matrix_2)

            J[:, i] = vector_2 / self.step
            print(datetime.now() - now)

        return J

    def drop_bad_singular_values(self, J: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        SVD decomposition of Jacobian. Some singular values can be dropped out.

        :param J: Jacobian
        :return: SVD decomposed Jacobian
        """
        svd = np.linalg.svd(np.matmul(J.T, J), full_matrices=False)
        u, sv, v = svd
        print("Singulars: ", sv)

        # plt.plot(sv, marker='o', markersize=5, label="Singular values of J")
        # plt.legend()
        # plt.show()

        sv = np.linalg.pinv(np.diag(sv))
        for i in range(len(sv)):
            if sv[i, i] > np.inf:
                sv[i, i] = 0

        return u, sv, v

    def calculate_parameters_delta(self, J: np.ndarray, u: np.ndarray, sv: np.ndarray, v: np.ndarray,
                                   vector_1: np.ndarray) -> np.ndarray:
        """
        Get final parameters deltas for elements to fit erroneous response matrix.

        :param J: Jacobian
        :param u: SVD of Jacobian
        :param sv: SVD of Jacobian
        :param v: SVD of Jacobian
        :param vector_1: interim residual vector
        :return: parameters deltas
        """
        J_new = np.matmul(np.matmul(v.T, sv), u.T)
        # delta = -np.matmul(np.linalg.pinv(np.matmul(J.T,J)),J.T).dot(vector_1)
        delta = -np.matmul(J_new, J.T).dot(vector_1)
        # delta = np.where(np.abs(delta)>1e-3, 0.5*delta, delta) # To remove twiss instabillity
        print("delta", delta)

        return delta

    def optimize_orbit(self) -> np.ndarray:
        """
        Calculate necessary parameters changes to make orbit correction.

        :return: parameters deltas for element alignments
        """
        experimental_responses = self.structure.collect_responses_on_quadrupole_kicks(self.structure.structure,
                                                                                      self.elements_to_vary,
                                                                                      gradient_step=1e-5)

        quadrupole_alignments = {}
        for i, quadrupole in enumerate(self.elements_to_vary):
            response_matrix, response_matrix_at_beginning = self.structure.calculate_response_matrix_on_quadrupole_variation(
                self.structure.structure,
                quadrupole,
                1e-5,
                i)
            delta = np.linalg.pinv(response_matrix).dot(experimental_responses[quadrupole[0]])
            quadrupole_alignments[quadrupole[0]] = delta
            print('Quadrupole alignment:', str(i) + '/' + str(self.elements_number))
            print('asdd', response_matrix_at_beginning)
            x = response_matrix_at_beginning.dot(delta)
            print("final", x)

        quadrupole_alignments = pd.DataFrame(quadrupole_alignments, index=['dx', 'dy'])

        return quadrupole_alignments


class LevenbergMarquardt(GaussNewton):
    def __init__(self, structure: Structure, step: float, iteration: int = 2, coefficient_lambda: float = 0.001):
        super().__init__(structure, step, iteration)
        self.coefficient_lambda = coefficient_lambda

    def drop_bad_singular_values(self, J: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        SVD decomposition of Jacobian. Some singular values can be dropped out.

        :param J: Jacobian
        :return: SVD decomposed Jacobian
        """
        svd = np.linalg.svd(np.matmul(J.T, J) + self.coefficient_lambda * np.diag(np.diag(np.matmul(J.T, J))),
                            full_matrices=False)
        u, sv, v = svd
        print("Singulars: ", sv)

        # plt.plot(sv, marker='o', markersize=5, label="Modified Jacobian")
        # svd_1 = np.linalg.svd(np.matmul(J.T, J), full_matrices=False)
        # _, sv_1, _ = svd_1
        # plt.plot(sv_1, color='r', marker='o', markersize=5, label="Not modified Jacobian")
        # plt.legend()
        # print("Singulars delta", sv - sv_1)
        # plt.show()

        sv = np.linalg.pinv(np.diag(sv))
        for i in range(len(sv)):
            if sv[i, i] > np.inf:
                sv[i, i] = 0

        return u, sv, v


class GaussNewtonConstrained(GaussNewton):
    def __init__(self, structure: Structure, step: float, iteration: int = 1):
        super().__init__(structure, step, iteration)
        self.weights = np.ones(self.elements_number) * 1e-3
        self.weights = np.random.rand(self.elements_number) * 1e-4

    def optimize_lattice(self) -> np.ndarray:
        """
        Calculate necessary parameters changes to make structure correction.

        :return: parameters deltas which fit erroneous structure
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
        initial_vector, initial_residual = self._get_residual(bad_response_matrix, model_response_matrix,
                                                              self.weights * delta)

        count = 1
        while count <= self.iteration:
            model_response_matrix_1 = self.structure.calculate_response_matrix(self.structure.structure,
                                                                               self.elements_to_vary,
                                                                               accumulative_param_additive,
                                                                               self.correctors)
            vector_1, _ = self._get_residual(bad_response_matrix, model_response_matrix_1,
                                             self.weights * delta)

            J = self.calculate_jacobian(accumulative_param_additive, model_response_matrix_1)
            u, sv, v = self.drop_bad_singular_values(J, delta)
            delta = self.calculate_parameters_delta(J_modified, u, sv, v, vector_1)
            accumulative_param_additive += delta
            count += 1

            fitted_model_response_matrix = self.structure.calculate_response_matrix(self.structure.structure,
                                                                                    self.elements_to_vary,
                                                                                    accumulative_param_additive,
                                                                                    self.correctors)
            final_vector, final_residual = self._get_residual(bad_response_matrix, fitted_model_response_matrix,
                                                              self.weights * delta)

            print("Initial parameters: ", self.initial_parameters)
            print("Bad parameters: ", self.bad_initial_parameters)
            print("Final parameters: ", list(self.initial_parameters + accumulative_param_additive))
            print("Initial residual: ", initial_residual)
            print("Final residual: ", final_residual)

        return accumulative_param_additive

    def create_modified_jacobian(self, J: np.ndarray, delta_qrad: np.ndarray) -> np.ndarray:
        """
        Function to get modified Jacobian.

        :param J: Jacobian
        :param delta_qrad: parameters deltas from last iteration
        :return: modified Jacobian
        """
        W = np.diag(self.weights * delta_qrad)
        return np.concatenate((J, W), axis=0)

    def drop_bad_singular_values(self, J: np.ndarray, delta_grad: np.ndarray) -> Tuple[
        np.ndarray, np.ndarray, np.ndarray]:
        """
        SVD decomposition of Jacobian. Some singular values can be dropped out.

        :param J: Jacobian
        :return: SVD decomposed Jacobian
        """
        J_modified = self.create_modified_jacobian(J, delta_qrad)
        svd = np.linalg.svd(np.matmul(J_modified.T, J_modified), full_matrices=False)
        u, sv, v = svd
        print("Singulars: ", sv)

        plt.plot(sv, marker='o', markersize=5, label="Modified Jacobian")
        svd_1 = np.linalg.svd(np.matmul(J.T, J), full_matrices=False)
        _, sv_1, _ = svd_1
        plt.plot(sv_1, color='r', marker='o', markersize=5, label="Not modified Jacobian")
        plt.legend()
        print("Singulars delta", sv - sv_1)
        plt.show()

        sv = np.linalg.pinv(np.diag(sv))
        for i in range(len(sv)):
            if sv[i, i] > np.inf:
                sv[i, i] = 0

        return u, sv, v
