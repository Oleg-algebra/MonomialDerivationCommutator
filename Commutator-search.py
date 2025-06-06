import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.odr import polynomial
from scipy.linalg import solve_banded,lstsq

from sympy import Matrix, solve_linear_system,symbols

class Monomial:
    def __init__(self, n,coeff: float,powers:list):
        self.n = n
        self.coeff = coeff
        if len(powers)!=n:
            raise AssertionError(f'powers must be of length {n}')
        self.powers = powers
    def getCoeff(self):
        return self.coeff
    def getPowers(self):
        return self.powers

    def __str__(self):
        monom = f"{self.coeff}"
        for i in range(len(self.powers)):
            monom += f" x_{i+1}^{self.powers[i]}"
        return monom

    def copy(self):
        return Monomial(self.n,self.coeff,self.powers.copy())


class Polynomial:
    def __init__(self, coefficients: np.ndarray,id):

        self.derivative_i_n = 0
        self.derivative_j_n = 0
        self.id = id
        self.coefficients = coefficients

    def copy(self):
        return Polynomial(self.coefficients.copy(),self.id)

    def __str__(self):
        return (f"Poly {self.id}, shape: {self.coefficients.shape}, der_i_n: {self.derivative_i_n}, der_j_n: {self.derivative_j_n}\n"
                f" coeffs: {self.coefficients}\n")

class MonomialCommutator:

    def __init__(self, n,monomials: list, N1, N2):
        for monomial in monomials:
            if len(monomial.getPowers()) != n:
                raise AssertionError(f'powers must be of length {n}')
        self.monomials = monomials
        self.n = n
        self.SOLE = None

        self.N1 = N1
        self.N2 = N2




    def generateCommutator(self) -> tuple[Polynomial,Polynomial]:

        powers1 = self.monomials[0].getPowers()
        powers2 = self.monomials[1].getPowers()
        N = abs(powers1[0] - powers2[0]) * self.N1 + self.N2
        M = abs(powers1[1] - powers2[1]) * self.N1 + self.N2
        coeffients = np.ones((N,M))
        return Polynomial(coeffients,0),Polynomial(coeffients,1)

    def x_derivative(self,polynomial: Polynomial) -> Polynomial:

        coefficients = polynomial.coefficients
        derivative_x = np.zeros(coefficients.shape)
        for i in range(1,coefficients.shape[0]):
            derivative_x[i-1] = coefficients[i] * i
        poly_derivative = Polynomial(derivative_x,polynomial.id)
        poly_derivative.derivative_i_n += 1
        return poly_derivative

    def y_derivative(self,polynomial:Polynomial) -> Polynomial:

        coefficients = polynomial.coefficients
        derivative_y = np.zeros(coefficients.shape)
        for i in range(1,coefficients.shape[1]):
            derivative_y[:,i-1] = coefficients[:,i] * i
        poly_derivative = Polynomial(derivative_y,polynomial.id)
        poly_derivative.derivative_j_n += 1

        return poly_derivative

    def monomial_derivative(self,monomial:Monomial,variable_n: int):
        if monomial.getPowers()[variable_n-1] == 0:
            derivative = monomial.copy()
            derivative.coeff = 0
        else:
            derivative = monomial.copy()
            derivative.coeff = monomial.coeff*monomial.getPowers()[variable_n-1]
            derivative.powers[variable_n-1] -= 1
        return derivative

    def create_BIG_matrix(self, polynomials: list[Polynomial],
                          monomial_derivatives: list[Monomial],
                          polynomial_derivatives: list[Polynomial],
                          ):

        N = max(self.monomials[0].getPowers()[0], self.monomials[1].getPowers()[0])
        M = max(self.monomials[0].getPowers()[1], self.monomials[1].getPowers()[1])

        start_i = []
        start_j = []
        end_i = []
        end_j = []

        BIG_Matrix = []
        for i in range(len(polynomials)):
            new_polynomial_coeff = np.zeros(shape=(polynomials[i].coefficients.shape[0] + N, polynomials[i].coefficients.shape[1] + M))

            T_monom = polynomials[i].coefficients * monomial_derivatives[i].coeff
            k1 = monomial_derivatives[i].powers[0]
            m1 = monomial_derivatives[i].powers[1]
            end11 = k1 + polynomials[i].coefficients.shape[0]
            end12 = m1 + polynomials[i].coefficients.shape[1]

            new_polynomial_coeff[k1:end11, m1:end12] = T_monom
            new_poly = Polynomial(new_polynomial_coeff,polynomials[i].id)
            BIG_Matrix.append(new_poly)

            start_i.append(k1)
            start_j.append(m1)
            end_i.append(end11)
            end_j.append(end12)

        for i in range(len(polynomials)):
            new_polynomial_coeff = np.zeros(shape=(polynomial_derivatives[i].coefficients.shape[0] + N, polynomial_derivatives[i].coefficients.shape[1] + M))

            D_poly = self.monomials[i].coeff * polynomial_derivatives[i].coefficients * (-1)

            # D_poly = Polynomial(D_poly,polynomial_derivatives[i].id)

            k3 = self.monomials[i].powers[0]
            m3 = self.monomials[i].powers[1]
            end31 = k3 + polynomial_derivatives[i].coefficients.shape[0]
            end32 = m3 + polynomial_derivatives[i].coefficients.shape[1]

            start_i.append(k3)
            start_j.append(m3)
            end_i.append(end31)
            end_j.append(end32)

            new_polynomial_coeff[k3:end31, m3:end32] = D_poly
            new_poly = Polynomial(new_polynomial_coeff, polynomial_derivatives[i].id)
            BIG_Matrix.append(new_poly)


        K1 = min(start_i)
        M1 = min(start_j)

        endK1 = max(end_i)
        endM1 = max(end_j)

        shifts = [[start_i[i],start_j[i]] for i in range(len(start_i))]

        return BIG_Matrix,[K1,endK1,M1,endM1],shifts

    def create_SOLE(self,polyPowers: list,
                    boundaries1: list,
                    boundaries2: list):
        cols = polyPowers[0]*polyPowers[1]*self.n
        rows = (boundaries1[1] - boundaries1[0] + 1)*(boundaries1[3] - boundaries1[2] + 1)
        rows += (boundaries2[1] - boundaries2[0] + 1)*(boundaries2[3] - boundaries2[2] + 1)

        return np.zeros((rows,cols))


    def fill_SOLE(self,BIG_Matrices: list[list[Polynomial]],
                  shifts_list: list,
                  boundaries_list: list,
                  polyPowers: list,
                  SOLE: np.ndarray):

        row_n = 0
        for m in range(len(BIG_Matrices)):
            matrices = BIG_Matrices[m]
            shifts = shifts_list[m]
            boundaries = boundaries_list[m]

            for i in range(boundaries[0],boundaries[1]):
                for j in range(boundaries[2],boundaries[3]):
                    row = np.zeros((1,polyPowers[0]*polyPowers[1]*2))
                    for k in range(len(matrices)):
                        a_ij = matrices[k].coefficients[i,j]

                        if a_ij == 0:
                            continue

                        ind_i = i - shifts[k][0] + matrices[k].derivative_i_n
                        ind_j = j - shifts[k][1] + matrices[k].derivative_j_n

                        position = matrices[k].id * (polyPowers[0]*polyPowers[1]) + ind_i * polyPowers[1] + ind_j

                        row[0,position] += a_ij

                    SOLE[row_n,:] = row
                    row_n += 1

        return SOLE

    def find_nonzero(self,SOLE,i):
        for j in range(i+1,SOLE.shape[0]):
            if SOLE[j,i] == 0:
                continue
            row_j = SOLE[j,:]
            row_i = SOLE[i,:]
            SOLE[i,:] = row_j
            SOLE[j,:] = row_i


    def direct_Gauss(self,SOLE: np.ndarray):

        for j in range(SOLE.shape[1]-1):
            row_j = SOLE[j,:]
            for i in range(j+1,SOLE.shape[0]):
                if SOLE[i,j] == 0:
                    continue
                if SOLE[j,j] == 0:
                    self.find_nonzero(SOLE,j)
                    if SOLE[j,j] == 0:
                        continue
                row_i = SOLE[i,:]
                lcm_ab = np.lcm(SOLE[j,j],SOLE[i,j])

                multiple1 = lcm_ab / SOLE[i,j]
                number = SOLE[j,j]
                multiple2 = lcm_ab / SOLE[j,j]
                m_row_i = row_i * multiple1
                m_row_j = row_j * multiple2
                new_row_i = m_row_i - m_row_j
                SOLE[i,:] = new_row_i

        return SOLE

    def backward_Gauss(self,SOLE: np.ndarray):

        for j in range(SOLE.shape[1]-1,0,-1):
            row_j = SOLE[j, :]
            for i in range(j - 1, -1,-1):
                if SOLE[i, j] == 0:
                    continue
                row_i = SOLE[i, :]
                lcm_ab = np.lcm(SOLE[j, j], SOLE[i, j])

                multiple1 = lcm_ab / SOLE[i, j]
                multiple2 = lcm_ab / SOLE[j, j]
                m_row_i = row_i * multiple1
                m_row_j = row_j * multiple2
                new_row_i = m_row_i - m_row_j
                SOLE[i, :] = new_row_i

        return SOLE


    def solve_SOLE(self,SOLE: np.ndarray):

        pass

    def commutator_search(self):
        [P,Q] = self.generateCommutator()
        # print("P shape:",P.coefficients.shape)
        # print("Q shape:",Q.coefficients.shape)
        P_derivatives = [self.x_derivative(P),self.y_derivative(P)]
        Q_derivatives = [self.x_derivative(Q),self.y_derivative(Q)]


        monomial1_derivatives = [self.monomial_derivative(self.monomials[0],1),
                                 self.monomial_derivative(self.monomials[0],2)]

        monomial2_derivatives = [self.monomial_derivative(self.monomials[1], 1),
                                 self.monomial_derivative(self.monomials[1], 2)]

        BIG_Matrix1,boundaries1,shifts1= self.create_BIG_matrix(
            [P,Q], monomial1_derivatives, P_derivatives,
        )
        BIG_Matrix2,boundaries2,shifts2 = self.create_BIG_matrix(
            [P,Q], monomial2_derivatives, Q_derivatives,
        )
        # print(shifts1)
        # print(boundaries1)
        # for poly in BIG_Matrix1:
        #     print(poly)

        # print(shifts2)
        # print(boundaries2)
        # print(BIG_Matrix2)
        SOLE = self.create_SOLE([P.coefficients.shape[0],P.coefficients.shape[1]],boundaries1,boundaries2)
        # print("SOLE shape:",SOLE.shape)
        SOLE = self.fill_SOLE([BIG_Matrix1,BIG_Matrix2],[shifts1,shifts2],
                              [boundaries1,boundaries2],P.coefficients.shape,SOLE)

        # print(SOLE)

        b = np.zeros((SOLE.shape[0],1))
        # print(b.shape)

        system_matrix = np.concatenate((SOLE,b),axis=1)
        # print(system_matrix.shape)
        A = Matrix(system_matrix)
        vars_a = []
        vars_b = []
        for i in range(P.coefficients.shape[0]):
            for j in range(P.coefficients.shape[1]):
                vars_a.append(symbols(f"a{i}_{j}"))
                vars_b.append(symbols(f"b{i}_{j}"))

        vars = np.concat((vars_a,vars_b))
        # print(vars)
        x = solve_linear_system(A,*vars)
        print(x)
        print(len(x.keys()))
        print(P.coefficients.shape)
        # print(A.shape)
        # print(len(vars))


def __str__(self):

        return "{} d/dx + {} d/dy".format(self.monomials[0],self.monomials[1])



if __name__ == "__main__":
    component1 = Monomial(2,-8,[6,7])
    component2 = Monomial(2,6,[5,8])
    moComm = MonomialCommutator(2,[component1,component2],1,2)

    moComm.commutator_search()
