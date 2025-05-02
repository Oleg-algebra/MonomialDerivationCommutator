import sympy
from numpy.polynomial.polynomial import Polynomial
from scipy.differentiate import derivative
from scipy.odr import polynomial
from sympy import Matrix, solve_linear_system,symbols, diff, simplify, expand, collect, Expr, solve, Poly, solve_linear, solve_undetermined_coeffs
import numpy as np

class Monomial:
    def __init__(self, n,coeff: float,powers:list):
        self.n = n
        self.coeff = coeff
        self.vars = []
        self.powers = powers
        if len(powers)!=n:
            raise AssertionError(f'powers must be of length {n}')
        self.monomial_symbolic = coeff
        for i in range(n):
            var = symbols(f"x_{i}")
            self.vars.append(var)
            self.monomial_symbolic *= var ** powers[i]
    def getCoeff(self):
        return self.coeff
    def getPowers(self):
        return self.powers

    def __str__(self):

        return self.monomial_symbolic



class Polynomial:
    def __init__(self, coefficients: np.ndarray = None,n_var: int = None,poly_symbols = None,vars = None):

        self.polynomial_symbolic = 0
        self.coefficients = coefficients
        if poly_symbols is not None:
            self.polynomial_symbolic = poly_symbols
        self.vars = []
        if vars is not None:
            self.vars = vars
        if coefficients is not None:
            for i in range(coefficients.shape[0]):
                for j in range(coefficients.shape[1]):
                    var = symbols(f"a{i}_{j}")
                    self.vars.append(var)
                    monomial = Monomial(n_var,coefficients[i,j],[i,j])
                    self.polynomial_symbolic += monomial.monomial_symbolic






class Derivation:
    def __init__(self, polynomials: list[Polynomial], variables):
        self.polynomials = polynomials
        self.variables = variables

    def take_derivative(self,expression: Expr ) -> Expr:
        res = 0
        for i in range(len(self.variables)):
            res += self.polynomials[i].polynomial_symbolic * diff(expression, self.variables[i])
        res = expand(res)
        return res


class Commutator:

    def __init__(self, derivation: Derivation,powers: list,K):
        self.derivation = derivation
        self.powers = powers
        self.unknown_coeffients = {}
        self.K = K

    def generateCommutator(self) -> Derivation:
        variables = self.derivation.variables
        N = max(abs(self.powers[0]-self.powers[1]),abs(self.powers[2]-self.powers[3])) + self.K
        print(N)
        Matrices = []
        symb = ["a","b"]

        for k in range(len(variables)):
            sym = symb[k]
            matrix = []
            for i in range(N):
                row = []
                for j in range(N):
                    coef = symbols(f'{sym}{i}_{j}')
                    row.append(coef)
                    self.unknown_coeffients[coef] = (i,j)
                matrix.append(row)
            Matrices.append(Matrix(matrix))
        polynomials = []
        for m in Matrices:
            polynomials.append(Polynomial(m,len(variables)))

        der_unknown = Derivation(polynomials,variables)

        return der_unknown

    def searchCommutator(self):
        unknown_derivation = self.generateCommutator()

        # for poly in self.derivation.polynomials:
        #     print(poly.polynomial_symbolic)
        #
        for poly in unknown_derivation.polynomials:
            print(f'unknown polynomial: {poly.polynomial_symbolic}')

        print(self.unknown_coeffients.keys())

        derivatives1 = []
        for poly in unknown_derivation.polynomials:
            derivatives1.append(self.derivation.take_derivative(poly.polynomial_symbolic))

        derivatives2 = []
        for poly in self.derivation.polynomials:
            derivatives2.append(unknown_derivation.take_derivative(poly.polynomial_symbolic))

        polys = []
        for i in range(len(derivatives1)):
            polys.append(derivatives1[i]-derivatives2[i])



        equations = []
        variables = unknown_derivation.variables
        for poly in polys:
            p = Poly(poly,variables)
            for term in p.terms():
                equations.append(term[1])
        print(f'equations: {equations}')
        res = solve(equations,list(self.unknown_coeffients.keys()))
        # res = solve_linear(equations,symbols=self.unknown_coeffients)
        res = solve_undetermined_coeffs(polys[0],coeffs=self.unknown_coeffients)
        res1 = solve_undetermined_coeffs(polys[1],coeffs=self.unknown_coeffients)
        return res, res1










if __name__ == "__main__":
    powers1 = [1,2]
    powers2 = [1,2]
    monomail1 = Monomial(2,-1,powers1)
    monomail2 = Monomial(2,1,powers2)
    polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
    polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)

    der = Derivation([polynomial1,polynomial2],polynomial2.vars)
    commutator = Commutator(der,[*powers1,*powers2],1)
    res, res1 = commutator.searchCommutator()
    print(res)
    print(res1)
    for coef in res.keys():
        print(coef," = ", res[coef])




