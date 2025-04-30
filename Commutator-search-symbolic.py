import sympy
from numpy.polynomial.polynomial import Polynomial
from scipy.differentiate import derivative
from scipy.odr import polynomial
from sympy import Matrix, solve_linear_system,symbols, diff, simplify, expand, collect, Expr
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

    def __init__(self, derivation: Derivation,powers: list):
        self.derivation = derivation
        self.powers = powers

    def generateCommutator(self) -> Derivation:
        variables = self.derivation.variables
        N = max(abs(self.powers[0]-self.powers[1]),abs(self.powers[2]-self.powers[3])) + 1
        print(N)
        Matrices = []
        symb = ["a","b"]

        for k in range(len(variables)):
            sym = symb[k]
            matrix = []
            for i in range(N):
                row = []
                for j in range(N):
                    row.append(f'{sym}{i}_{j}')
                matrix.append(row)
            Matrices.append(Matrix(matrix))
        polynomials = []
        for m in Matrices:
            polynomials.append(Polynomial(m,len(variables)))

        der_unknown = Derivation(polynomials,variables)

        return der_unknown

    def searchCommutator(self):
        unknown_derivation = self.generateCommutator()
        for poly in self.derivation.polynomials:
            print(poly.polynomial_symbolic)

        for poly in unknown_derivation.polynomials:
            print(poly.polynomial_symbolic)

        derivatives1 = []
        for poly in unknown_derivation.polynomials:
            derivatives1.append(self.derivation.take_derivative(poly.polynomial_symbolic))

        derivatives2 = []
        for poly in self.derivation.polynomials:
            derivatives2.append(unknown_derivation.take_derivative(poly.polynomial_symbolic))

        equations = []
        for i in range(len(derivatives1)):
            equations.append(derivatives1[i]-derivatives2[i])


        new_equations = []
        for eq in equations:
            new_eq = eq.subs(symbols("x_0"),symbols("x"))
            new_eq = new_eq.subs(symbols("x_1"),symbols("y"))
            new_equations.append(new_eq)

        x = symbols("x")
        y = symbols("y")
        for eq in new_equations:
            print(eq)

        #TODO: write code to collect coefficients by x^i * y^j
        print(new_equations[0].coeff(x**2*y**2))







if __name__ == "__main__":
    powers1 = [1,2]
    powers2 = [1,2]
    monomail1 = Monomial(2,1,powers1)
    monomail2 = Monomial(2,1,powers2)
    polynomial1 = Polynomial(poly_symbols=monomail1.monomial_symbolic,vars=monomail1.vars)
    polynomial2 = Polynomial(poly_symbols=monomail2.monomial_symbolic,vars=monomail2.vars)

    der = Derivation([polynomial1,polynomial2],polynomial2.vars)
    commutator = Commutator(der,[*powers1,*powers2])
    commutator.searchCommutator()



