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
# Performs random tests to help create a hypothesis about commutator nature
import os
import time

from CommutatorSearchSymbolic import *

def do_tests(
        tests_number = 10,
        powers = ([1,1],[1,1]),
        coefficients = (1,1),
        variables_number = 2,
        K = 2 ,      # additional bias for number of unknown coefficients for polynomials
        fileName = os.path.basename(__file__).split(".")[0]

):


    proportionalCounter = 0         # counts solutions that proportional to D
    unproportionalCounter = 0       #  counts solutions that are not proportional to D
    zeroDerivationCounter = 0

    results = {}
    givenDerivationKEY = "GIVEN_DERIVATION"
    commutatorKEY = "COMMUTATOR"
    proportionalKEY = "proportionalCounter"
    unproportionalKEY = "unproportionalCounter"
    isProportionalKEY = "isProportional"
    isZeroDerivationKEY = "IsZeroDerivation"
    matrixDimension  = "matrixDimension"

    s = 0

    for i in range(tests_number):
        start = time.time()

        print("Testing " + str(i+1) + "/" + str(tests_number))
        result = {}

        powers1 = powers[0]
        powers2 = powers[1]

        alpha = coefficients[0]
        beta = coefficients[1]


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

        zeroCounter = 0
        for i in range(len(res.polynomials)):
            commutatorPolynomials.append(res.polynomials[i].polynomial_symbolic)
            if res.polynomials[i].polynomial_symbolic.equals(0):
                zeroCounter += 1
        if zeroCounter == variables_number:
            zeroDerivationCounter += 1
            result[isZeroDerivationKEY] = True

        else:
            result[isZeroDerivationKEY] = False



        result[commutatorKEY] = commutatorPolynomials
        result[isProportionalKEY] = isProportional


        if isProportional:
            if zeroCounter != variables_number:
                proportionalCounter +=1
        else:
            unproportionalCounter+=1

        results[(*powers1,*powers2)] = result


        end = time.time()
        s+=(end-start)



    results[proportionalKEY] = proportionalCounter
    results[unproportionalKEY] = unproportionalCounter
    results[zeroDerivationCounter] = zeroDerivationCounter


    print(f'proportional: {proportionalCounter}')
    print(f'unproportional: {unproportionalCounter}')
    print(f'zeroDerivaions: {zeroDerivationCounter}')
    print("Total time: ", s)
    print("avarage time: ", s / tests_number)

    file = open(fileName+"_log.txt", "w")

    file.write("Report of testing\n")
    file.write("======================General information===================\n")
    file.write(f"Number of tests: {tests_number}\n")
    file.write(f"proportional: {proportionalCounter}\n")
    file.write(f"unproportional: {unproportionalCounter}\n")
    file.write(f"zeroDerivaions: {zeroDerivationCounter}\n")
    file.write(f"total time: {s}\n")
    file.write(f"avarage time: {s / tests_number}\n")
    file.write("======================Special cases=========================\n")

    file.write("=====================Proportional derivations==================\n")
    for param, res in results.items():
        if type(res).__name__ == "dict":
            if res[isProportionalKEY] and not res[isZeroDerivationKEY]:
                file.write(f"{param}: {res}\n")

    file.write("=======================Zero derivations=========================\n")
    for param, res in results.items():
        if type(res).__name__ == "dict":
            if res[isZeroDerivationKEY]:
                file.write(f"{param}: {res}\n")

    file.write("=======================Unproportional derivation=========================\n")
    for param, res in results.items():
        if type(res).__name__ == "dict":
            if not res[isProportionalKEY]:
                file.write(f"{param}: {res}\n")


    file.write("==================END of REPORT=======================\n")
    file.close()