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


from CommutatorSearchSymbolic import *

tests_number = 1000

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
isZeroDerivationKEY = "IsZeroDerivation"
matrixDimension  = "matrixDimension"

s = 0
K = 2
counter = 0

while counter < tests_number:
    start = time.time()

    print("Testing " + str(counter+1) + "/" + str(tests_number))
    result = {}
    l = np.random.randint(0, max_power)
    k = l + 1
    n = np.random.randint(0, max_power)
    m = n + 1

    alpha = np.random.randint(min_coeff, max_coeff)
    beta = np.random.randint(min_coeff, max_coeff)
    # beta = alpha

    a = np.random.randint(min_coeff, max_coeff)
    alpha = -a * m
    beta = a * k

    powers1 = [k,n]
    powers2 = [l, m]

    if alpha ** 2 + beta **2 == 0:
        alpha = np.random.randint(min_coeff, max_coeff)
        beta = np.random.randint(min_coeff, max_coeff)

    monomail1 = Monomial(variables_number,alpha,powers1)
    monomail2 = Monomial(variables_number,beta,powers2)
    polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
    polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)



    result[givenDerivationKEY] = [polynomial1.polynomial_symbolic, polynomial2.polynomial_symbolic]


    der = Derivation([polynomial1,polynomial2],monomail1.vars)

    commutator = Commutator(der,[*powers1,*powers2],K)
    print(f'Matrix dimension: {commutator.unknown_derivation.polynomials[0].coefficients.shape}')


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
        result[isZeroDerivationKEY] = True
    else:
        result[isZeroDerivationKEY] = False

    result[commutatorKEY] = commutatorPolynomials
    result[isProportionalKEY] = isProportional

    if (k,n,l,m) not in results.keys():
        results[(k,n,l,m)] = result
        counter += 1
    else:
        continue
    print(len(results.keys()))

    end = time.time()
    s+=(end-start)

for res in results.values():
    if res[isZeroDerivationKEY]:
        zeroDerivationCounter +=1
    if res[isProportionalKEY] and not res[isZeroDerivationKEY]:
        proportionalCounter += 1
    if not res[isProportionalKEY] :
        unproportionalCounter += 1

results[proportionalKEY] = proportionalCounter
results[unproportionalKEY] = unproportionalCounter
results[zeroDerivationCounter] = zeroDerivationCounter

print(len(results.keys()))

print(f'proportional: {proportionalCounter}')
print(f'unproportional: {unproportionalCounter}')
print(f'zeroDerivaions: {zeroDerivationCounter}')
print("Total time: ", s)
print("avarage time: ", s / tests_number)

fileName = os.path.basename(__file__).split(".")[0]
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
count = 1
for param, res in results.items():
    if type(res).__name__ == "dict":
        if res[isProportionalKEY] and not res[isZeroDerivationKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1

file.write("=======================Zero derivations=========================\n")
count = 1
for param, res in results.items():
    if type(res).__name__ == "dict":
        if res[isZeroDerivationKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1

file.write("=======================Unproportional derivation=========================\n")
count = 1
for param, res in results.items():
    if type(res).__name__ == "dict":
        if not res[isProportionalKEY]:
            file.write(f"{count}):  {param}: {res}\n")
            count += 1


file.write("==================END of REPORT=======================\n")
file.close()

