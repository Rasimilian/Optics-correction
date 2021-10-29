import numpy as np
import pandas as pd

from cpymad.madx import Madx


def read_elements_from_file(elements_file):
    """
    Read elements list from file.

    :param str elements_file: file name with elements
    :return: list of elements, amount of elements
    """
    elements = open(elements_file, 'r')

    count = 1
    elements_list = []
    for element in elements:
        elements_list.append(element.strip())
        count += 1
    elements.close()

    return elements_list, count


def describe_correctors(structure: str, elements_file: str):
    """
    Find predefined elements in structure and write their parameters.

    :param list structure_in_lines: opened structure in lines
    :return: list of elements parameters, list of the last parameters
    """
    elements, count = read_elements_from_file(elements_file)

    elements_description = []
    elements_parameters = []
    elements = [i.lower() for i in elements]

    madx = Madx(stdout=False)
    madx.option(echo=False, warn=False, info=False, twiss_print=False)
    madx.call(file=structure)
    madx.input('beam,particle=electron,energy=1.8;')

    madx_elements = madx.elements
    for i, element in enumerate(madx_elements):
        if element in elements:
            elements_description.append((element, i, madx_elements[i].kick))
            elements_parameters.append(madx_elements[i].kick)
    madx.quit()

    return elements_description, elements_parameters


def describe_elements(structure: str, elements_file: str):
    """
    Find predefined elements in structure and write their parameters.

    :param list structure_in_lines: opened structure in lines
    :return: list of elements parameters, list of the last parameters
    """
    elements, count = read_elements_from_file(elements_file)

    elements_description = []
    elements_parameters = []
    elements = [i.lower() for i in elements]

    madx = Madx(stdout=False)
    madx.option(echo=False, warn=False, info=False, twiss_print=False)
    madx.call(file=structure)
    madx.input('beam,particle=electron,energy=1.8;')

    madx_elements = madx.elements
    for i, element in enumerate(madx_elements):
        if element in elements:
            elements_description.append((element,i,madx_elements[i].k1))
            elements_parameters.append(madx_elements[i].k1)
    madx.quit()

    return elements_description, elements_parameters


def describe_elements1(structure_in_lines):
    """
    Find predefined elements in structure and write their parameters.

    :param list structure_in_lines: opened structure in lines
    :return: list of elements parameters, list of the last parameters
    """
    elements, count = read_elements_from_file("madx\elements\quads.txt")
    # print(elements)

    elements_description = []
    elements_parameters = []

    for element in elements:
        element = element.split()[0]
        for line in structure_in_lines:
            if line.startswith(element):
                line = line.replace(',', '').replace(';', '').split()
                elements_description.append(line)
                elements_parameters.append(float(line[-1]))

    elements_parameters = np.array(elements_parameters)

    # print(elements_description)
    return elements_description, elements_parameters


def read_BPMs(BPMs_data):
    """
    Read data from BPMs.

    :param str BPMs_data: file name with BPMs data
    :return: float X, Y beam orbits
    """
    data = pd.read_csv(BPMs_data, sep='\t')
    assert data.shape != 54

    X = data.iloc[:,1:3].to_numpy()
    Y = data.iloc[:,1:4:2].to_numpy()

    return X, Y


def add_errors_to_structure_file(structure='madx\structures\VEPP4M_full1.txt', elements='madx\elements\elems.txt'):
    elements_description, elements_parameters = describe_elements(structure, elements)

    madx = Madx(stdout=False)
    madx.option(echo=False, warn=False, info=False, twiss_print=False)
    madx.call(file=structure)
    madx.input('beam,particle=electron,energy=1.8;')
    madx.input('use,sequence=RING;')

    for element in elements_description:
        print(madx.elements[element[1]].k1)
        madx.elements[element[1]].k1 = 1.01*element[2]
        print(madx.elements[element[1]].k1)

    madx.input('save, sequence=RING, file="VEPP4M_full1_all_errors1112.txt"')
    madx.quit()

add_errors_to_structure_file()

