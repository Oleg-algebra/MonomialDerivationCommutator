from CommutatorSearchSymbolic import *
from math import sqrt

variables_number = 2
isSearchNonZero = True
isZeroDerivation = False
max_K = 25

k = 9
n= 14
l = 9
m = 14

alpha = -33
beta = -21

powers1 = [k, n]
powers2 = [l, m]

monomail1 = Monomial(2,alpha,powers1)
monomail2 = Monomial(2,beta,powers2)
polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)
print(polynomial1.variables_polynom)

der = Derivation([polynomial1,polynomial2], polynomial1.variables_polynom)
K = 1
commutator = Commutator(der,[*powers1,*powers2],K)

print("========Given derivation=======")
for i in range(len(der.polynomials)):
    print(f'poly {i}: {der.polynomials[i].polynomial_symbolic}')

print("==========Unknown derivation=======")
while True:
    print(f"K = {K}")
    commutatorPolynomials = []
    commutator = Commutator(der, [*powers1, *powers2], K)
    print(f"Matrices size: {sqrt(len(commutator.unknown_coeffients) / 2)}")

    res, isProportional = commutator.get_commutator()


    zeroCounter = 0
    for i in range(len(res.polynomials)):
        commutatorPolynomials.append(res.polynomials[i].polynomial_symbolic)
        if res.polynomials[i].polynomial_symbolic.equals(0):
            zeroCounter += 1
    if zeroCounter == variables_number:
        if K < max_K and isSearchNonZero:
            print("---> zeroDerivation")
            print("="*100)
            K += 1
            continue
        isZeroDerivation = True
    else:
        # if K<max_K and isProportional:
        #     print("--->nonZeroDerivation and Proportional")
        #     for i in range(len(res.polynomials)):
        #         print(f'poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}')
        #         print("="*100)
        #     K += 1
        #     continue
        isZeroDerivation = False
    break

for i in range(len(res.polynomials)):
    print(f'poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}' )

print(f"proportional: {isProportional}")
print(f"is Solution correct: {commutator.isSolution(derivation1=der,derivation2=res)}")
print(f'isZeroCommutator: {isZeroDerivation}')
