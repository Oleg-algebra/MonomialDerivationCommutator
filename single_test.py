from CommutatorSearchSymbolic import *

k = 0
n= 9
l = 0
m = 9

alpha = 23
beta = -26

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
res, isProportional = commutator.get_commutator()
for i in range(len(res.polynomials)):
    print(f'poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}' )

print(f"proportional: {isProportional}")
print(f"is Solution correct: {commutator.isSolution(derivation1=der,derivation2=res)}")
