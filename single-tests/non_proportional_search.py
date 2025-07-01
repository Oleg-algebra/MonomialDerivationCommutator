from CommutatorSearchSymbolic import *
from math import sqrt
import argparse
import os
from time import time
from tqdm import tqdm

variables_number = 2
isSearchNonZero = True
isZeroDerivation = False

startTime = time()

parser = argparse.ArgumentParser(description="A simple script with command-line arguments.")
parser.add_argument("--N", help="power")
parser.add_argument("--K", help="max bias power")

args = parser.parse_args()
# N = int(args.N)
# max_K = int(args.K)
N = 1
max_K = 3
k = 0
n= N
l = N
m = 0

k = 2
n = 6
l = 2
m = 6

alpha = 8
beta = 0

powers1 = [k, n]
powers2 = [l, m]

monomail1 = Monomial(2,alpha,powers1)
monomail2 = Monomial(2,beta,powers2)
polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)


der = Derivation([polynomial1,polynomial2], polynomial1.variables_polynom)
K = 2
commutator = Commutator(der,[*powers1,*powers2],K)
file = open(f"test_{N}_output.txt","w")
file.write("Test results\n")
file.write("="*100 + "\n")
file.write(f"Given derivation: {[polynomial1.polynomial_symbolic,polynomial2.polynomial_symbolic]}\n")

with tqdm(total=max_K,desc=f"Case N = {N}",position = 0,leave=False,disable=False) as pbar:
    pbar.update(2)
    while K < max_K:
        # print(f"K = {K}")
        commutatorPolynomials = []
        commutator = Commutator(der, [*powers1, *powers2], K)
        # print(f"Matrices size: {commutator.unknown_derivation.polynomials[0].coefficients.shape}")

        res, isProportional = commutator.get_commutator()

        zeroCounter = 0
        for i in range(len(res.polynomials)):
            commutatorPolynomials.append(simplify(res.polynomials[i].polynomial_symbolic))
            if res.polynomials[i].polynomial_symbolic.equals(0):
                zeroCounter += 1

        if zeroCounter == 2:
            isZeroDerivation = True
        else:
            isZeroDerivation = False

        file.write(f"K = {K} --- isProportional: {isProportional} --- COMMUTATOR: {commutatorPolynomials}\n")

        if isProportional:
            # print("Proportional ---> increase K by 1")
            K += 1
        else:
            break

        pbar.update(1)

endTime = time()
totalTime = endTime - startTime
file.close()

print(f"Tests for N = {N} finished")
print(f"Total time: {totalTime}")
print(f"Variables: {polynomial1.variables_polynom}")
print("========Given derivation=======")
for i in range(len(der.polynomials)):
    print(f'poly {i}: {der.polynomials[i].polynomial_symbolic}')

print("==========Unknown derivation=======")

print(f"proportional: {isProportional}")
print(f"is Solution correct: {commutator.isSolution(derivation1=der,derivation2=res)}")
print("FINAL COMMUTATOR")
for i in range(len(res.polynomials)):
    print(f'poly {i}: {simplify(res.polynomials[i].polynomial_symbolic)}' )
print("="*100)
reportFile = open(f"report.txt","a")
if isProportional:
    reportFile.write(f"{N} 1\n")
else:
    reportFile.write(f"{N} 0\n")
reportFile.close()

