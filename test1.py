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
import time

import numpy as np
from CommutatorSearchSymbolic import *

tests_number = 500
true_counter = 0        # counts solutions that proportional to D
false_counter = 0       #  counts solutions that are not proportional to D

min_coeff = -10
max_coeff = 10

max_power = 10

variables_number = 2
proportionalCounter = 0
unproportionalCounter = 0
zeroDerivationCounter = 0

results = {}
givenDerivationKEY = "GIVEN_DERIVATION"
commutatorKEY = "COMMUTATOR"
proportionalKEY = "proportionalCounter"
unproportionalKEY = "unproportionalCounter"
isProportionalKEY = "isProportional"
zeroDerivaionsKEY = "zeroDerivaions"

s = 0

for i in range(tests_number):
    start = time.time()

    print("Testing " + str(i+1) + "/" + str(tests_number))
    result = {}
    l = np.random.randint(0, max_power)
    k = l + 1
    n = np.random.randint(0, max_power)
    m = n + 1

    # l = m = k = n = 1
    # alpha = -1
    # beta = 1

    powers1 = [k,n]
    powers2 = [l, m]

    alpha = np.random.randint(min_coeff, max_coeff)
    beta = np.random.randint(min_coeff, max_coeff)
    # beta = alpha

    a = np.random.randint(min_coeff, max_coeff)
    alpha = -a*m
    beta = a * k

    if alpha ** 2 + beta **2 == 0:
        print("alpha == beta == 0")
        print("----> SKIP ---->")
        continue

    monomail1 = Monomial(variables_number,alpha,powers1)
    monomail2 = Monomial(variables_number,beta,powers2)
    polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
    polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)



    result[givenDerivationKEY] = [polynomial1.polynomial_symbolic, polynomial2.polynomial_symbolic]


    der = Derivation([polynomial1,polynomial2],polynomial2.vars)
    K = 2
    commutator = Commutator(der,[*powers1,*powers2],K)





    res, isProportional = commutator.get_commutator()
    commutatorPolynomials = []
    if isProportional and alpha != -beta:
        print("-->degrees: ", (k, n, l, m))
        print(f"-->alpha: {alpha}, beta: {beta}")
        print("========Given derivation=======")
        for i in range(len(der.polynomials)):
            print(f'-->poly {i}: {der.polynomials[i].polynomial_symbolic}')

        print("==========Unknown derivation=======")
        zeroCounter = 0
        for i in range(len(res.polynomials)):
            print(f'-->poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}' )
            commutatorPolynomials.append(res.polynomials[i].polynomial_symbolic)
            if res.polynomials[i].polynomial_symbolic.equals(0):
                zeroCounter += 1
            print(res.polynomials[i].polynomial_symbolic.equals(0))
        if zeroCounter == variables_number:
            zeroDerivationCounter += 1
    result[commutatorKEY] = commutatorPolynomials
    result[isProportionalKEY] = isProportional

    print(f"proportional: {isProportional}")

    if isProportional:
        proportionalCounter +=1
    else:
        unproportionalCounter+=1

    results[(k,n,l,m)] = result


    end = time.time()
    print("Time elapsed: " + str(end-start))
    s+=(end-start)
    print("=" * 200)



results[proportionalKEY] = proportionalCounter
results[unproportionalKEY] = unproportionalCounter
results[zeroDerivationCounter] = zeroDerivationCounter

# for key,values in results.items():
#     print(key)
#     print(values)
#     print("================================")
print(f'proportional: {proportionalCounter}')
print(f'unproportional: {unproportionalCounter}')
print(f'zeroDerivaions: {zeroDerivationCounter}')
print("Total time: ", s)
print("avarage time: ", s / tests_number)