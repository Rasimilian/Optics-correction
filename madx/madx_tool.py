import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from typing import List, Tuple
from cpymad.madx import Madx

from element_parser.data_parser import read_elements_from_file, describe_elements


class Structure():
    def __init__(self,
                 structure_file: str = "madx\structures\VEPP4M_full1.txt",
                 bad_structure_file: str = "madx\structures\VEPP4M_full1_all_grads_errors.txt"):
        """
        Initialize class Structure.

        :param str structure_file: file name with structure
        :param str bad_structure_file: file name with bad structure
        """
        self.structure = structure_file
        self.structure_in_lines = open(self.structure).readlines()
        self.bad_structure = bad_structure_file
        self.bad_structure_in_lines = open(self.bad_structure).readlines()

        self.BPMs_number = 54
        # self.elements_number = len(describe_elements(self.structure_in_lines)[0])

        self.twiss_table_4D = self.calculate_structure_4D(self.structure)
        self.twiss_table_6D = self.calculate_structure_6D(self.structure)
        self.bad_twiss_table_4D = self.calculate_structure_4D(self.bad_structure)
        self.bad_twiss_table_6D = self.calculate_structure_6D(self.bad_structure)

    def calculate_structure_4D(self, structure: str, initial_imperfections=None) -> Madx.twiss:
        """
        Calculate TWISS table for 4D beam motion.

        :param str structure: file name with structure
        :param initial_imperfections:
        :return: float TWISS table
        """
        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')
        madx.input('use,sequence=RING;')

        # TODO add errors
        # self.errors_table = self.add_errors_to_structure(madx, 'quadrupole', value=0.000001, initial_imperfections=initial_imperfections)

        # madx.input('select,flag=twiss,class=monitor;')

        madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_BPMs;')
        twiss_table = madx.table.twiss_in_BPMs
        # madx.quit()

        return twiss_table

    def calculate_structure_6D(self, structure: str, initial_imperfections=None) -> Madx.twiss:
        """
        Calculate TWISS table for 6D beam motion.

        :param str structure: file name with structure
        :param initial_imperfections:
        :return: float TWISS table
        """
        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')
        madx.input('use,sequence=RING;')

        # TODO add errors
        # self.errors_table = self.add_errors_to_structure(madx, 'quadrupole', value=0.000001, initial_imperfections=initial_imperfections)

        # madx.input('select,flag=twiss,class=monitor;')
        madx.input('ptc_create_universe;ptc_create_layout,model=2,method=2,nst=1;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=6,nst=10,exact=true;')

        madx.ptc_twiss(icase=6, no=1, center_magnets=True, closed_orbit=True, table='twiss', file="madx\\log_file.txt")
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_BPMs;')

        twiss_table = madx.table.twiss_in_BPMs
        # madx.quit()

        return twiss_table

    def make_kick_by_corrector(self,
                               corrector: Tuple[str, int, float],
                               madx: Madx,
                               kick_step: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Change elements parameters in structure.

        :param str structure: file name with structure
        :param list structure_in_lines: opened structure in lines
        :param int parameter_number: number of varying parameter
        :param float variation_step: step to vary elements parameters
        :param float accumulative_param_additive: to accumulate parameters changes after iterations
        :return: float TWISS table
        """
        # now = datetime.now()
        madx.elements[corrector[1]].kick = corrector[2] + kick_step

        # madx.input('select,flag=twiss,class=monitor;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=2,nst=1;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=6,nst=10,exact=true;')
        # madx.ptc_twiss(icase=6,no=1,center_magnets=True,closed_orbit=True,table='twiss', file="madx\\log_file.txt")

        madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
        # madx.input('twiss,sequence=RING, centre=True, table=twiss, file=madx\\log_file.txt;')
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_bpms;')

        twiss_table = madx.table.twiss_in_bpms.x, madx.table.twiss_in_bpms.y
        # print(twiss_table[0])
        madx.elements[corrector[1]].kick = corrector[2]
        # print(datetime.now()-now)
        return twiss_table

    def make_kick_by_quadrupole(self,
                                quadrupole: Tuple[str, int, float],
                                madx: Madx,
                                gradient_step: float) -> Tuple[np.ndarray, np.ndarray]:
        madx.elements[quadrupole[1]].k1 = quadrupole[2] + gradient_step

        madx.input('select,flag=twiss,class=monitor;')
        madx.input('ptc_create_universe;ptc_create_layout,model=2,method=2,nst=1;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=6,nst=10,exact=true;')

        madx.ptc_twiss(icase=6, no=1, center_magnets=True, closed_orbit=True, table='twiss',
                       file="madx\\log_file.txt", )

        # madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
        # madx.input('twiss,sequence=RING, centre=True, table=twiss, file=madx\\log_file.txt;')
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_bpms;')

        twiss_table = madx.table.twiss_in_bpms.x, madx.table.twiss_in_bpms.y
        madx.elements[quadrupole[1]].k1 = quadrupole[2]

        return twiss_table

    def change_structure_for_correction(self,
                                        structure: str,
                                        elements_to_vary: List[Tuple[str, int, float]],
                                        accumulative_param_additive: np.ndarray) -> Madx.twiss:
        """
        Change elements parameters in structure for correction. All locations are involved.

        :param str structure: file name with structure
        :param list structure_in_lines: opened structure in lines
        :param int parameter_number: number of varying parameter
        :param float variation_step: step to vary elements parameters
        :param float accumulative_param_additive: to accumulate parameters changes after iterations
        :return: float TWISS table
        """
        number_of_elements = len(elements_to_vary)
        accumulative_param_additive_ = accumulative_param_additive

        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')

        for num, element in enumerate(elements_to_vary):
            madx.elements[element[1]].k1 = element[2] + accumulative_param_additive_[num]

        madx.input('use,sequence=RING;')

        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=2,nst=1;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=6,nst=10,exact=true;')
        # madx.ptc_twiss(icase=6,no=1,center_magnets=True,closed_orbit=True,table='twiss', file="madx\\log_file.txt",)
        madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_BPMs;')

        twiss_table = madx.table.twiss_in_BPMs
        # madx.quit()

        return twiss_table

    def calculate_response_matrix(self,
                                  structure: str,
                                  elements_to_vary: List[Tuple[str, int, float]],
                                  accumulative_param_additive: np.ndarray,
                                  correctors: List[Tuple[str, int, float]],
                                  corrector_step: float = 1e-6,
                                  accumulative_alignment_additive=None,
                                  areErrorsNeeded=None) -> pd.DataFrame:
        """
        Calculate response matrix.

        :param str structure: file name with structure
        :param list structure_in_lines: opened structure in lines
        :param float variation_step: step to vary elements parameters
        :param float accumulative_param_additive: to accumulate parameters changes after iterations
        :param float accumulative_alignment_additive: to accumulate alignment errors changes after iterations
        :param bool areErrorsNeeded: whether to add errors
        :return: float response matrix
        """
        # now = datetime.now()
        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')
        madx.input('use,sequence=RING;')
        madx.input('select,flag=twiss,class=monitor;')

        accumulative_param_additive_ = accumulative_param_additive
        for num, element in enumerate(elements_to_vary):
            madx.elements[element[1]].k1 = element[2] + accumulative_param_additive_[num]

        frames = []
        twiss_in_BPMs = self.make_kick_by_corrector(correctors[0], madx, kick_step=0)
        for corrector in correctors:
            ## Jacobian = (f(x+dx)-f(x))/dx
            ## For f(x+dx)
            twiss_in_BPMs_1 = self.make_kick_by_corrector(corrector, madx, corrector_step)

            ## For f(x)
            # twiss_in_BPMs = self.measure_response(self.structure, self.structure_in_lines, element, variation_step[0],
            #                                     accumulative_param_additive[n])
            # print(len(twiss_in_BPMs))
            # print(len(twiss_in_BPMs_1))
            # assert len(twiss_in_BPMs_1) != len(twiss_in_BPMs)

            df_x = pd.DataFrame((twiss_in_BPMs_1[0] - twiss_in_BPMs[0]) / corrector_step, columns=[corrector[0]])
            df_y = pd.DataFrame((twiss_in_BPMs_1[1] - twiss_in_BPMs[1]) / corrector_step, columns=[corrector[0]])
            # df_dx = pd.DataFrame((twiss_in_BPMs_1.dx - twiss_in_BPMs.dx)/variation_step, columns = [element[0]])
            # df_dy = pd.DataFrame((twiss_in_BPMs_1.dy - twiss_in_BPMs.dy)/variation_step, columns = [element[0]])

            df = pd.concat([df_x, df_y], ignore_index=True)
            # df = pd.concat([df_x,df_y,df_dx,df_dy], ignore_index = True)
            frames.append(df)

        matrix = pd.concat(frames, axis=1)
        # print("Response Matrix:\n",matrix)
        # matrix.to_csv('madx//response_matrix_quad_x.txt',index=False,header=False,sep="\t")
        # print("resp:",datetime.now()-now)
        madx.quit()
        return matrix

    def calculate_response_matrix_on_quadrupole_variation(self,
                                                          structure: str,
                                                          quadrupole: Tuple[str, int, float],
                                                          gradient_step: float,
                                                          num: int,
                                                          alignment_step: float = 1e-6) -> Tuple[pd.DataFrame, pd.DataFrame]:
        alignments = ['dx', 'dy']

        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')
        madx.input('use,sequence=RING;')
        madx.input('eoption,add=False;')

        twiss_in_BPMs = self.calculate_twiss_in_bpms(madx, False)
        frames = []
        madx.input('select,flag=error,pattern=' + quadrupole[0] + ';')
        for alignment in alignments:
            madx.input('ealign,' + alignment + ':=' + str(alignment_step) + ';')
            twiss_in_BPMs_1 = self.calculate_twiss_in_bpms(madx, False)
            df_x = pd.DataFrame([(twiss_in_BPMs_1[0][num] - twiss_in_BPMs[0][num]) / alignment_step],
                                columns=[quadrupole[0] + '_' + alignment])
            df_y = pd.DataFrame([(twiss_in_BPMs_1[1][num] - twiss_in_BPMs[1][num]) / alignment_step],
                                columns=[quadrupole[0] + '_' + alignment])
            df = pd.concat([df_x, df_y], ignore_index=True)
            frames.append(df)
        madx.input('ealign,' + alignment + ':=' + str(0) + ';')
        matrix0 = pd.concat(frames, axis=1)

        twiss_in_BPMs = self.calculate_twiss_in_bpms(madx, True)
        frames = []
        madx.input('select,flag=error,pattern=' + quadrupole[0] + ';')
        for alignment in alignments:
            madx.input('ealign,' + alignment + ':=' + str(alignment_step) + ';')
            twiss_in_BPMs_1 = self.calculate_twiss_in_bpms(madx, True)
            df_x = pd.DataFrame((twiss_in_BPMs_1[0] - twiss_in_BPMs[0]) / alignment_step,
                                columns=[quadrupole[0] + '_' + alignment])
            df_y = pd.DataFrame((twiss_in_BPMs_1[1] - twiss_in_BPMs[1]) / alignment_step,
                                columns=[quadrupole[0] + '_' + alignment])
            df = pd.concat([df_x, df_y], ignore_index=True)
            frames.append(df)
        madx.input('ealign,' + alignment + ':=' + str(0) + ';')

        madx.input('select,flag=errors_table,class=quadrupole;')
        madx.input('etable,table=errors_table;')
        print('asd', madx.table.errors_table.dx, madx.table.errors_table.dy)
        matrix = pd.concat(frames, axis=1)

        frames = []
        madx.elements[quadrupole[1]].k1 = quadrupole[2] + gradient_step
        twiss_in_BPMs = self.calculate_twiss_in_bpms(madx, True)
        for alignment in alignments:
            madx.input('ealign,' + alignment + ':=' + str(alignment_step) + ';')
            twiss_in_BPMs_1 = self.calculate_twiss_in_bpms(madx, True)
            df_x = pd.DataFrame((twiss_in_BPMs_1[0] - twiss_in_BPMs[0]) / alignment_step,
                                columns=[quadrupole[0] + '_' + alignment])
            df_y = pd.DataFrame((twiss_in_BPMs_1[1] - twiss_in_BPMs[1]) / alignment_step,
                                columns=[quadrupole[0] + '_' + alignment])
            df = pd.concat([df_x, df_y], ignore_index=True)
            frames.append(df)
        matrix1 = pd.concat(frames, axis=1)
        madx.input('select,flag=errors_table,class=quadrupole;')
        madx.input('etable,table=errors_table;')
        print(madx.table.errors_table.dx, madx.table.errors_table.dy)

        madx.quit()
        return matrix1 - matrix, matrix0

    def calculate_twiss_in_bpms(self, madx: Madx, at_centre: bool) -> Tuple[np.ndarray, np.ndarray]:
        madx.input('select,flag=twiss,clear;')
        if at_centre == True:
            madx.input('select,flag=twiss,class=monitor;')
        else:
            madx.input('select,flag=twiss,pattern=^mq;')
        # madx.input('ptc_create_universe;ptc_create_layout,model=2,method=2,nst=1;')
        #
        # madx.ptc_twiss(icase=6,no=1,center_magnets=True,closed_orbit=True,table='twiss', file="madx\\log_file.txt",)

        # madx.twiss(sequence='RING', centre=True, table='twiss', file="madx\\log_file.txt")
        madx.input('twiss,sequence=RING, centre=False, table=twiss, file=madx\\log_file.txt;')
        madx.input('readtable,file="madx\\log_file.txt",table=twiss_in_bpms;')

        return madx.table.twiss_in_bpms.x, madx.table.twiss_in_bpms.y

    def collect_responses_on_quadrupole_kicks(self, structure, elements_to_vary, gradient_step, alignment_error=1e-4):
        madx = Madx(stdout=False)
        madx.option(echo=False, warn=False, info=False, twiss_print=False)
        madx.call(file=structure)
        madx.input('beam,particle=electron,energy=1.8;')
        madx.input('use,sequence=RING;')

        madx.input('select,flag=error,pattern=' + 'SIL1' + ';')
        madx.input('select,flag=error,pattern=' + 'SIL2' + ';')
        madx.input('select,flag=error,pattern=' + 'STL1' + ';')
        madx.input('ealign,dx:=' + str(alignment_error) + ',dy:=' + str(alignment_error) + ';')

        madx.input('select,flag=errors_table,class=quadrupole;')
        madx.input('etable,table=errors_table;')

        twiss_in_BPMs = self.calculate_twiss_in_bpms(madx, True)
        frames = []
        for i, quadrupole in enumerate(elements_to_vary):
            twiss_in_BPMs_1 = self.make_kick_by_quadrupole(quadrupole, madx, gradient_step)
            df_x = pd.DataFrame((twiss_in_BPMs_1[0] - twiss_in_BPMs[0]), columns=[quadrupole[0]])
            df_y = pd.DataFrame((twiss_in_BPMs_1[1] - twiss_in_BPMs[1]), columns=[quadrupole[0]])
            df = pd.concat([df_x, df_y], ignore_index=True)
            frames.append(df)
        matrix = pd.concat(frames, axis=1)

        madx.quit()
        return matrix

    def add_imperfections_to_structure(self,
                                       elements_with_errors,
                                       accumulative_alignment_additive,
                                       error_magnitude=1e-6):
        seed = np.random.randint(0, 999999999)
        self.madx.input('eoption,seed=' + str(seed) + ',add=True;')

        for element, dx, dy in zip(elements_with_errors, accumulative_alignment_additive,
                                   accumulative_alignment_additive):
            madx.input('select,flag=error,pattern=' + element + ';')
            madx.input('ealign,dx:=' + str(error_magnitude) + '*gauss(),dy:=' + str(error_magnitude) + '*gauss();')
            madx.input('select,flag=error,clear;')

        madx.input('esave,file="madx\machine_imperfections";')
        # madx.input('select,flag=myerrortable, class=quadrupole;')
        madx.input('etable,table=errors_table;')

        errors_table = madx.table.errors_table


class Imperfection(Structure):
    def __init__(self, madx, types_of_errors, file_with_elements_to_spoil):
        """
        Initialize class Imperfection.

        :param object madx: MAD-X interpreter
        :param list types_of_errors: list of errors types
        :param str file_with_elements_to_spoil: file with listed elements to spoil
        """
        self.madx = madx
        self.elements_list, self.elements_number = read_elements_from_file(file_with_elements_to_spoil)

        if not isinstance(types_of_errors, list): self.types_of_errors = list[types_of_errors]

    def add_errors(self, error_amplitude=0.000001, element_type='quadrupole'):
        """
        Add errors to magnetic elements.

        :param float error_amplitude: magnitude of errors
        :param str element_type: type of elements to add errors to
        :return: float table with created errors
        """
        seed = np.random.randint(0, 999999999)
        self.madx.input('eoption,seed=' + str(seed) + ',add=True;')

        for element in elements:
            # TODO check the line below
            element = element.strip()
            madx.input('select,flag=error,class=' + element_type + ',pattern=' + element + ';')

        # TODO add types_of_errors choice
        madx.input('ealign,dx:=' + str(error_amplitude) + '*gauss(),dy:=' + str(error_amplitude) + '*gauss();')
        madx.input('esave,file="madx\machine_imperfections";')
        madx.input('select,flag=myerrortable, class=quadrupole;')
        madx.input('etable,table=errors_table;')

        errors_table = madx.table.errors_table

        return errors_table
