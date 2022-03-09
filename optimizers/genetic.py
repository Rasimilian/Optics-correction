# import random
# import math
# import matplotlib.pyplot as plt
#
#
# class GA(object):
#     # Инициализировать популяцию Создание популяции из людей размера_хромосомы размера_длины_хромосомы
#
#     def __init__(self, population_size, chromosome_length, max_value, pc, pm):
#
#         self.population_size = population_size
#         self.choromosome_length = chromosome_length
#         # self.population=[[]]
#         self.max_value = max_value
#         self.pc = pc
#         self.pm = pm
#         # self.fitness_value=[]
#
#     def species_origin(self):
#         population = [[]]
#         for i in range(self.population_size):
#
#             temporary = []
#             # Регистр хромосом
#             for j in range(self.choromosome_length):
#                 temporary.append(random.randint(0, 1))
#             # Произвольно сгенерировать хромосому, состоящую из двоичных чисел
#
#             population.append(temporary)
#             # Добавить хромосомы в популяцию
#         return population[1:]
#         # Возвращаем популяцию, популяция - это двумерный массив, индивидуум и хромосома - двумерные
#
#     # От двоичного к десятичному
#     # Кодирование ввода: популяция, длина хромосомы Процесс кодирования - это процесс преобразования многомерной функции в унарную функцию
#     def translation(self, population):
#
#         temporary = []
#         for i in range(len(population)):
#             total = 0
#             for j in range(self.choromosome_length):
#                 total += population[i][j] * (math.pow(2, j))
#             # Начиная с первого гена, каждый человек возводится в степень до 2, а затем суммируется
#             # Например: 0101 преобразовано в десятичное: 1 * 20 + 0 * 21 + 1 * 22 + 0 * 23 = 1 + 0 + 4 + 0 = 5
#             temporary.append(total)
#         # Завершено кодирование хромосомы, которое преобразуется двоичным числом в десятичное.
#         return temporary
#
#     # Возвращает десятичное число всех особей в популяции после кодирования
#
#     # from protein to function,according to its functoin value
#
#     # a protein realize its function according its structure
#     # Целевая функция эквивалентна среде для фильтрации хромосом, здесь 2 * sin (x) + math.cos (x)
#     def function(self, population):
#         temporary = []
#         function1 = []
#         temporary = self.translation(population)
#         for i in range(len(temporary)):
#             x = temporary[i] * self.max_value / (math.pow(2, self.choromosome_length) - 10)
#             y =
#             function1.append(-(x-1)**2)
#             # (x**2 + y - 11)**2 + (x + y**2 -7)**2
#
#         # Здесь sin (x) используется как целевая функция
#         return function1
#
#     # Определить фитнес
#     def fitness(self, function1):
#
#         fitness_value = []
#
#         num = len(function1)
#
#         for i in range(num):
#
#             if (function1[i] < 0):
#                 temporary = function1[i]
#             else:
#                 temporary = 0.0
#             # Если фитнес меньше 0, он устанавливается на 0
#
#             fitness_value.append(temporary)
#         # Добавить фитнес в список
#
#         return fitness_value
#
#     # Рассчитать фитнес и
#
#     def sum(self, fitness_value):
#         total = 0
#
#         for i in range(len(fitness_value)):
#             total += fitness_value[i]
#         return total
#
#     # Рассчитать фитнес Fibe и список
#     def cumsum(self, fitness1):
#         for i in range(len(fitness1) - 2, -1, -1):
#             # range(start,stop,[step])
#             # Обратный отсчет
#             total = 0
#             j = 0
#
#             while (j <= i):
#                 total += fitness1[j]
#                 j += 1
#
#             fitness1[i] = total
#             fitness1[len(fitness1) - 1] = 1
#
#
#     def selection(self, population, fitness_value):
#         new_fitness = []
#         # Регистр единой формулы
#         total_fitness = self.sum(fitness_value)
#         # Подсчитайте все фитнес
#         for i in range(len(fitness_value)):
#             new_fitness.append(fitness_value[i] / total_fitness)
#         # Упорядочить физическую форму всех людей
#         self.cumsum(new_fitness)
#         #
#         ms = []
#         # Выжившее население
#         population_length = pop_len = len(population)
#         # Найти длину популяции
#         # Определить, какие из них могут выжить, по случайным числам
#
#         for i in range(pop_len):
#             ms.append(random.random())
#         # Генерация случайного значения численности населения
#         # ms.sort()
#         # Сортировка выжившего населения
#         fitin = 0
#         newin = 0
#         new_population = new_pop = population
#
#         # Путь рулетки
#         while newin < pop_len:
#             if (ms[newin] < new_fitness[fitin]):
#                 new_pop[newin] = population[fitin]
#                 newin += 1
#             else:
#                 fitin += 1
#         population = new_pop
#
#
#     def crossover(self, population):
#         # pc - это порог вероятности, выберите одноточечный кроссовер или многоточечный кроссовер для создания новых кроссоверов, что здесь бесполезно
#         pop_len = len(population)
#
#         for i in range(pop_len - 1):
#
#             if (random.random() < self.pc):
#                 cpoint = random.randint(0, len(population[0]))
#                 # Произвольно генерировать единственную точку пересечения в пределах количества популяций
#                 temporary1 = []
#                 temporary2 = []
#
#                 temporary1.extend(population[i][0:cpoint])
#                 temporary1.extend(population[i + 1][cpoint:len(population[i])])
#                 # Использовать tmporary1 как временную память для временного хранения первых генов от 0 до cpoint в i-й хромосоме,
#                 # Затем добавьте последнюю точку c в i + 1-й хромосоме к количеству генов в i-й хромосоме и добавьте ее в конец временного 2
#
#                 temporary2.extend(population[i + 1][0:cpoint])
#                 temporary2.extend(population[i][cpoint:len(population[i])])
#                 # Используйте tmporary2 как временную память для временного хранения первых генов от 0 до cpoint в i + 1-й хромосоме,
#                 # Затем добавьте последнюю точку c в i-й хромосоме к количеству генов в i-й хромосоме и добавьте ее в конец временного 2
#                 population[i] = temporary1
#                 population[i + 1] = temporary2
#         # Рекомбинация / кроссовер гена i-й хромосомы и i + 1-й хромосомы завершена
#
#
#     def mutation(self, population):
#         # pm - порог вероятности
#         px = len(population)
#         # Найдите количество всех популяций / особей в популяции
#         py = len(population[0])
#         # Количество хромосом / отдельных генов
#         for i in range(px):
#             if (random.random() < self.pm):
#                 mpoint = random.randint(0, py - 1)
#                 #
#                 if (population[i][mpoint] == 1):
#                     # Произвести одноточечную случайную мутацию генов mpoint до 0 или 1
#                     population[i][mpoint] = 0
#                 else:
#                     population[i][mpoint] = 1
#
#
#     # transform the binary to decimalism
#     # Преобразуйте каждую хромосому в десятичное значение max_value, а затем отфильтруйте лишнее значение
#     def b2d(self, best_individual):
#         total = 0
#         b = len(best_individual)
#         for i in range(b):
#             total = total + best_individual[i] * math.pow(2, i)
#
#         total = total * self.max_value / (math.pow(2, self.choromosome_length) - 1)
#         return total
#
#
#     # Найдите лучший фитнес и индивидуальный
#
#     def best(self, population, fitness_value):
#         px = len(population)
#         bestindividual = []
#         bestfitness = fitness_value[0]
#         # print(fitness_value)
#
#         for i in range(1, px):
#             # Цикл, чтобы найти максимальную физическую форму, лучший человек - лучший человек
#             if (fitness_value[i] > bestfitness):
#                 bestfitness = fitness_value[i]
#                 bestindividual = population[i]
#
#         return [bestindividual, bestfitness]
#
#
#     def plot(self, results):
#         X = []
#         Y = []
#
#         for i in range(500):
#             X.append(i)
#             Y.append(results[i][0])
#
#         plt.plot(X, Y)
#         plt.show()
#
#
#     def main(self):
#         results = [[]]
#         fitness_value = []
#         fitmean = []
#
#         population = pop = self.species_origin()
#
#         for i in range(500):
#             function_value = self.function(population)
#             # print('fit funtion_value:',function_value)
#             fitness_value = self.fitness(function_value)
#             # print('fitness_value:',fitness_value)
#
#             best_individual, best_fitness = self.best(population, fitness_value)
#             results.append([best_fitness, self.b2d(best_individual)])
#             # Сохраните лучшие индивидуальные и лучшие фитнес-функции и преобразуйте лучших в десятичные числа, функция фитнеса
#             self.selection(population, fitness_value)
#             self.crossover(population)
#             self.mutation(population)
#         results = results[1:]
#         results.sort()
#         self.plot(results)
#
#
# if __name__ == '__main__':
#     population_size = 1000
#     max_value = 10
#     chromosome_length = 20
#     pc = 0.6
#     pm = 0.01
#     ga = GA(population_size, chromosome_length, max_value, pc, pm)
#     ga.main()


from deap import base, algorithms, creator, tools

import random
import matplotlib.pyplot as plt
import numpy as np


# LOW, UP = -5, 5
LOW, UP = -1e-4, 1e-4
ETA = 20
LENGTH_CHROM = 26
POPULATION_SIZE = 10
P_CROSSOVER = 0.9
P_MUTATION = 0.2
MAX_GENERATIONS = 50
HALL_OF_FAME_SIZE = 5

hof = tools.HallOfFame(HALL_OF_FAME_SIZE)

RANDOM_SEED = 42
random.seed(RANDOM_SEED)

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

from madx.madx_tool import Structure
from optimizers import GaussNewton
structure = Structure()
optimizer = GaussNewton(structure, "madx/correctors/correctors.txt", "madx/elements/quads.txt", 1e-6,
                        grad_step=1e-3)
accumulative_param_additive = np.zeros(optimizer.elements_number)
accumulative_alignment_additive = {'dx': np.zeros(optimizer.elements_number), 'dy': np.zeros(optimizer.elements_number)}

real_matrix = structure.calculate_response_matrix(optimizer.structure.bad_structure,
                                                  optimizer.bad_elements_to_vary,
                                                  accumulative_param_additive,
                                                  optimizer.bad_correctors,
                                                  optimizer.corrector_step,
                                                  optimizer.structure.quads_with_alignments,
                                                  accumulative_alignment_additive)

def randomPoint(a, b):
    return np.random.uniform(a, b, size=optimizer.elements_number)


toolbox = base.Toolbox()
toolbox.register('randomPoint', randomPoint, LOW, UP)
toolbox.register("individualCreator", tools.initIterate, creator.Individual, toolbox.randomPoint)
toolbox.register("populationCreator", tools.initRepeat, list, toolbox.individualCreator)

population = toolbox.populationCreator(n=POPULATION_SIZE)


def himmeblau(individual):
    print(individual)
    matrix = structure.calculate_response_matrix(optimizer.structure.structure,
                                                 optimizer.elements_to_vary,
                                                 individual,
                                                 optimizer.correctors,
                                                 optimizer.corrector_step,
                                                 optimizer.names,
                                                 accumulative_alignment_additive)
    _, res = optimizer._get_residual(real_matrix, matrix)
    print(res)
    return res,

# def himmeblau(individual):
#     x, y = individual
    # return (x**2 + y - 11)**2 + (x + y**2 -7)**2,
    # return -(y+47)*np.sin(np.sqrt(np.abs(x/2+y+47)))-x*np.sin(np.sqrt(np.abs(x-y-47))),

from deap import tools
from deap.algorithms import varAnd


def eaSimpleElitism(population, toolbox, cxpb, mutpb, ngen, stats=None,
             halloffame=None, verbose=__debug__, callback=None):
    """Перелеланный алгоритм eaSimple с элементом элитизма
    """

    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    hof_size = len(halloffame.items) if halloffame.items else 0

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population) - hof_size)

        # Vary the pool of individuals
        offspring = varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        offspring.extend(halloffame.items)

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        # Replace the current population by the offspring
        population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)

        if callback:
            callback[0](*callback[1])

    return population, logbook


toolbox.register('evaluate', himmeblau)
toolbox.register('select', tools.selTournament, tournsize=3)
# toolbox.register('mate', tools.cxSimulatedBinaryBounded, low=LOW, up=UP, eta=ETA)
# toolbox.register('mutate', tools.mutPolynomialBounded, low=LOW, up=UP, eta=ETA, indpb=1.0/LENGTH_CHROM)
toolbox.register('mate', tools.cxSimulatedBinary, eta=ETA)
toolbox.register('mutate', tools.mutFlipBit, indpb=1.0/LENGTH_CHROM)

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("min", np.min)
stats.register("avg", np.mean)

import time
def show(ax, xgrid, ygrid, f):
    ptMins = [[3.0, 2.0], [-2.805118, 3.131312], [-3.779310, -3.283186], [3.584458, -1.848126]]

    ax.clear()
    ax.contour(xgrid, ygrid, f)
    ax.scatter(*zip(*ptMins), marker='X', color='red', zorder=1)
    ax.scatter(*zip(*population), color='green', s=2, zorder=0)

    plt.draw()
    plt.gcf().canvas.flush_events()

    time.sleep(0.2)


x = np.arange(-5, 5, 0.1)
y = np.arange(-5, 5, 0.1)
xgrid, ygrid = np.meshgrid(x, y)

f_himmelblau = (xgrid**2 + ygrid - 11)**2 + (xgrid + ygrid**2 - 7)**2

# plt.ion()
# fig, ax = plt.subplots()
# fig.set_size_inches(5, 5)
#
# ax.set_xlim(LOW-3, UP+3)
# ax.set_ylim(LOW-3, UP+3)

#algelitism.eaSimpleElitism
#algorithms.eaSimple
population, logbook = eaSimpleElitism(population, toolbox,
                                        cxpb=P_CROSSOVER,
                                        mutpb=P_MUTATION,
                                        ngen=MAX_GENERATIONS,
                                        halloffame=hof,
                                        stats=stats,
                                        # callback=(show, (ax, xgrid, ygrid, f_himmelblau)),
                                        verbose=True)

maxFitnessValues, meanFitnessValues = logbook.select("min", "avg")

best = hof.items[0]
print(best)

# plt.ioff()
# plt.show()

# plt.plot(maxFitnessValues, color='red')
# plt.plot(meanFitnessValues, color='green')
# plt.xlabel('Поколение')
# plt.ylabel('Макс/средняя приспособленность')
# plt.title('Зависимость максимальной и средней приспособленности от поколения')
# plt.show()
