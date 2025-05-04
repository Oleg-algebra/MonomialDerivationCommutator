#
#
# D = alpha*x^k*y^n* d/dx + beta*x^l*y^m d/dy  -- given derivation
#
#T = P(x,y)d/dx + Q(x,y)d/dy  -- unknown derivation
#
#
#where P = sum_i,j = 0 ^ n a_ij*x^iy^j
#       Q = sum_i,j = 0 ^ n b_ij*x^iy^j
#     b_ij, a_ij  --- needs to be found
#
#
# case 1: l=k and n = m
import os
import time
from mpi4py import MPI
from tqdm import tqdm
import sys
from cases_functions import get_parameters


from CommutatorSearchSymbolic import *

sys.setrecursionlimit(10**6)

def isSolution(derivation1: Derivation,derivation2: Derivation) -> bool:
    polyDerivatives1 = []
    polyDerivatives2 = []

    for poly in derivation1.polynomials:
        der = derivation2.take_derivative(poly.polynomial_symbolic)
        polyDerivatives1.append(der)

    for poly in derivation2.polynomials:
        der = derivation1.take_derivative(poly.polynomial_symbolic)
        polyDerivatives2.append(der)

    for i in range(len(polyDerivatives1)):
        difference = polyDerivatives1[i] - polyDerivatives2[i]
        if not difference.equals(0):
            return False
    return True

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

case = 1

tests_number = 100

tests_number = tests_number // size + 2

coeff_limit = 50
min_coeff = -coeff_limit
max_coeff = coeff_limit

max_power = 10
min_power = 0

K = 2
max_K = 7

variables_number = 2
proportionalCounter = 0
unproportionalCounter = 0
zeroDerivationCounter = 0
correct_answers_counter = 0
false_answers_counter = 0

results = {}
givenDerivationKEY = "GIVEN_DERIVATION"
commutatorKEY = "COMMUTATOR"
isProportionalKEY = "isProportional"
isZeroDerivationKEY = "IsZeroDerivation"
matrixDimension  = "matrixDimension"
isSolutionCorrectKey = "isSolutionCorrect"

correctAnswersNumberKEY = "correctAnswersNumber"
falseAnswersNumberKey = "falseAnswersNumber"
proportionalKEY = "proportionalCounter"
unproportionalKEY = "unproportionalCounter"
zeroDerivaionsKEY = "zeroDerivaionsCounter"
time_exec_KEY = "time_elapsed"

s = 0

l = 0
k = 0
n = 0
m = 0

alpha = 0
beta = 0

counter = 0
# print(f'rank: {rank} started testing')
with tqdm(total=tests_number,desc=f"Rank: {rank}",position=rank,leave=False) as pbar:
    while counter < tests_number:

        start = time.time()

        # print(f"rank {rank} ---> Testing {str(counter+1) + "/" + str(tests_number)}" )
        result = {}

        #   Determine parameters
        #########################################################################
        if K == 2:
            # print("--->")
            l = np.random.randint(0, max_power)
            k = l
            n = np.random.randint(0, max_power)
            m = n

            alpha = np.random.randint(min_coeff, max_coeff)
            beta = np.random.randint(min_coeff, max_coeff)

            l,k,n,m,alpha,beta = get_parameters(case, min_power, max_power, min_coeff, max_coeff)
        else:
            # print("--> increased K by 1")
            pass
        ############################################################################

        powers1 = [k,n]
        powers2 = [l, m]

        if alpha ** 2 + beta **2 == 0:
            continue

        if (k,n,l,m,alpha,beta) in results.keys():
            continue

        monomail1 = Monomial(variables_number,alpha,powers1)
        monomail2 = Monomial(variables_number,beta,powers2)
        polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
        polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)

        result[givenDerivationKEY] = [polynomial1.polynomial_symbolic, polynomial2.polynomial_symbolic]

        der = Derivation([polynomial1,polynomial2],monomail1.vars)

        commutator = Commutator(der,[*powers1,*powers2],K)

        result[isZeroDerivationKEY] = False

        res, isProportional = commutator.get_commutator()
        result[matrixDimension] = commutator.unknown_derivation.polynomials[0].coefficients.shape
        commutatorPolynomials = []

        if isSolution(der,res):
            result[isSolutionCorrectKey] = True
        else:
            result[isSolutionCorrectKey] = False

        zeroCounter = 0
        for i in range(len(res.polynomials)):
            commutatorPolynomials.append(res.polynomials[i].polynomial_symbolic)
            if res.polynomials[i].polynomial_symbolic.equals(0):
                zeroCounter += 1
        if zeroCounter == variables_number:
            if K < max_K:
                K += 1
                continue
            result[isZeroDerivationKEY] = True

        else:

            result[isZeroDerivationKEY] = False

        result["K"] = K
        result[commutatorKEY] = commutatorPolynomials
        result[isProportionalKEY] = isProportional


        end = time.time()
        time_elapsed = end - start
        result[time_exec_KEY] = time_elapsed

        s+=time_elapsed
        results[(k, n, l, m, alpha, beta)] = result
        K = 2
        counter += 1
        pbar.update(1)


comm.Barrier()
results_container = comm.gather(results, root=0)
comm.Barrier()
if rank == 0:

    average_K = 0
    max_K = 2

    total_time = 0
    all_results = {}
    for dct in results_container:
        for key in dct.keys():
            all_results[key] = dct[key]

    for res in all_results.values():

        average_K += res["K"]
        if res["K"] > max_K:
            max_K = res["K"]

        if res[isZeroDerivationKEY]:
            zeroDerivationCounter += 1
        if res[isProportionalKEY] and not res[isZeroDerivationKEY]:
            proportionalCounter += 1
        if not res[isProportionalKEY]:
            unproportionalCounter += 1
        if res[isSolutionCorrectKey]:
            correct_answers_counter += 1
        else:
            false_answers_counter += 1
        total_time += res[time_exec_KEY]

    average_K = average_K / len(all_results.keys())
    average_time_per_process = total_time / size

    print("\n"+"="*100)
    print(f"Total number of different cases: {len(all_results.keys())}")
    print(f"rank: {rank} gathering results")
    # print(results_container)

    print(f'proportional: {proportionalCounter}')
    print(f'unproportional: {unproportionalCounter}')
    print(f'zeroDerivations: {zeroDerivationCounter}')
    print(f"Average K: {average_K}")
    print(f"Max K: {max_K}")
    print("Average time per process: ", average_time_per_process)
    print("average time per test: ", total_time / (tests_number*size))
    print("="*100)

    fileName = f"case_{case}"
    file = open(fileName+"_log.txt", "w")

    file.write("Report of testing\n")
    file.write("======================General information===================\n")
    file.write(f"Total number of different cases: {len(all_results.keys())}\n")
    file.write(f"proportional: {proportionalCounter}\n")
    file.write(f"unproportional: {unproportionalCounter}\n")
    file.write(f"zeroDerivations: {zeroDerivationCounter}\n")
    file.write(f"correct answers number: {correct_answers_counter}\n")
    file.write(f"false answers number: {false_answers_counter}\n")
    file.write(f"Average K: {average_K}")
    file.write(f"Max K: {max_K}")
    file.write(f"Average time per process: {average_time_per_process}")
    file.write(f"average time per test: {total_time / (tests_number*size)}\n")
    file.write("======================Special cases=========================\n")

    file.write("=====================Proportional derivations==================\n")
    count = 1
    for param, res in all_results.items():

        if res[isProportionalKEY] and not res[isZeroDerivationKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1

    file.write("=======================Zero derivations=========================\n")
    count = 1
    for param, res in all_results.items():
        if res[isZeroDerivationKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1

    file.write("=======================Unproportional derivation=========================\n")
    count = 1
    for param, res in all_results.items():
        if not res[isProportionalKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1

    file.write("=======================Correct answers=========================\n")
    count = 1
    for param, res in all_results.items():
        if res[isSolutionCorrectKey]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1


    file.write("=======================False answers=========================\n")
    count = 1
    for param, res in all_results.items():
        if not res[isSolutionCorrectKey]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1



    file.write("==================END of REPORT=======================\n")
    file.close()


