from numpy.polynomial.polynomial import Polynomial
from sympy import Matrix, solve_linear_system,symbols, diff, simplify, expand, Expr, solve, Poly, N, nsimplify
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
        self.variables_polynom = []
        if vars is not None:
            self.variables_polynom = vars
        if coefficients is not None:
            for i in range(coefficients.shape[0]):
                for j in range(coefficients.shape[1]):
                    monomial = Monomial(n_var,coefficients[i,j],[i,j])
                    if self.variables_polynom == []:
                        self.variables_polynom = monomial.vars
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
        self.unknown_coeffients = []
        self.K = K
        self.unknown_derivation = self.generateCommutator()
        self.searchCommutator = {
            "general" : self.generalSolver,
            "linear" : self.linearSolver
        }

    def generateCommutator(self) -> Derivation:
        variables = self.derivation.variables
        N = max(abs(self.powers[0]-self.powers[1]),abs(self.powers[2]-self.powers[3])) + self.K

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
                    self.unknown_coeffients.append(coef)
                matrix.append(row)
            Matrices.append(Matrix(matrix))
        polynomials = []
        for m in Matrices:
            polynomials.append(Polynomial(m,len(variables)))

        der_unknown = Derivation(polynomials,variables)

        return der_unknown

    def generalSolver(self):

        # for poly in self.unknown_derivation.polynomials:
        #     print(f'unknown polynomial: {poly.polynomial_symbolic}')

        # print(self.unknown_coeffients)

        derivatives1 = []
        for poly in self.unknown_derivation.polynomials:
            derivatives1.append(self.derivation.take_derivative(poly.polynomial_symbolic))

        derivatives2 = []
        for poly in self.derivation.polynomials:
            derivatives2.append(self.unknown_derivation.take_derivative(poly.polynomial_symbolic))

        polys = []
        for i in range(len(derivatives1)):
            polys.append(derivatives1[i]-derivatives2[i])



        equations = []
        variables = self.unknown_derivation.variables
        for poly in polys:
            p = Poly(poly,variables)
            for term in p.terms():
                equations.append(term[1])
        # print(f'equations: {equations}')
        res = solve(equations,self.unknown_coeffients)

        return res

    def linearSolver(self):


        derivatives1 = []
        for poly in self.unknown_derivation.polynomials:
            derivatives1.append(self.derivation.take_derivative(poly.polynomial_symbolic))

        derivatives2 = []
        for poly in self.derivation.polynomials:
            derivatives2.append(self.unknown_derivation.take_derivative(poly.polynomial_symbolic))

        polys = []
        for i in range(len(derivatives1)):
            polys.append(derivatives1[i] - derivatives2[i])

        equations = []
        variables = self.unknown_derivation.variables

        terms = []
        for poly in polys:
            p = Poly(poly, variables)
            for coeff in p.terms():
                p1  = Poly(coeff[1],self.unknown_coeffients)
                terms.append(p1)

        for term in terms:
            # print(term)
            row = np.zeros((1,len(self.unknown_coeffients)))
            for t in term.terms():
                coeffs = np.array(t[0])
                scalar = float(N(t[1], chop=True))
                row = row + coeffs * scalar
            equations.append(row[0])

        # print(f'equations: {equations}')
        matrix = np.array(equations)
        # print(matrix)
        matrix = np.concatenate((matrix,np.zeros((matrix.shape[0],1))),axis=1)
        matrix = Matrix(matrix)
        res = solve_linear_system(matrix, *self.unknown_coeffients)

        return res

    def get_commutator(self,solver = "general"):

        coefficients = self.searchCommutator[solver]()

        arbitrary_coefficients = []
        for coeff in self.unknown_coeffients:
            if coeff not in coefficients.keys():
                arbitrary_coefficients.append(coeff)


        for poly in self.unknown_derivation.polynomials:
            for coeff in coefficients.keys():
                new_symbolic_expr = poly.polynomial_symbolic.subs(coeff,coefficients[coeff])
                poly.polynomial_symbolic = new_symbolic_expr



        for poly in self.unknown_derivation.polynomials:
            for coeff in arbitrary_coefficients:
                number = np.random.randint(1,10)
                new_symbolic_expr = poly.polynomial_symbolic.subs(coeff,1)
                poly.polynomial_symbolic = nsimplify(new_symbolic_expr,rational=True)
                # poly.polynomial_symbolic = new_symbolic_expr

        is_proportional = self.is_proportional()
        return self.unknown_derivation,is_proportional


    def is_proportional(self):
        poly_unknown = self.unknown_derivation.polynomials
        poly_given = self.derivation.polynomials

        fractions = []



        for i in range(len(poly_unknown)):
            fraction = poly_unknown[i].polynomial_symbolic / poly_given[i].polynomial_symbolic
            fraction = simplify(fraction)
            fractions.append(fraction)
        # print(fractions)
        const = simplify(fractions[0]/fractions[0])
        for i in range(1,len(fractions)):
            check  = fractions[0].equals(fractions[i])
            if not check:
                return False

        return True


    def isSolution(self,derivation1: Derivation,derivation2: Derivation) -> bool:
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








