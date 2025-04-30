import numpy as np

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
    def __init__(self, coefficients):
        self.coefficients = coefficients

    def copy(self):
        return Polynomial(self.coefficients.copy())

class MonomialCommutator:

    def __init__(self, n,monomials: list):
        for monomial in monomials:
            if len(monomial.getPowers()) != n:
                raise AssertionError(f'powers must be of length {n}')
        self.monomials = monomials
        self.n = n
        self.SOLE = None




    def generateCommutator(self) -> tuple[Polynomial,Polynomial]:

        powers1 = self.monomials[0].getPowers()
        powers2 = self.monomials[1].getPowers()
        N = abs(powers1[0] - powers2[0]) + 1
        M = abs(powers1[1] - powers2[1]) + 1
        coeffients = np.ones((N,M))
        return Polynomial(coeffients),Polynomial(coeffients)

    def x_derivative(self,polynomial: Polynomial) -> Polynomial:

        coefficients = polynomial.coefficients
        derivative_x = np.zeros(coefficients.shape)
        for i in range(1,coefficients.shape[0]):
            derivative_x[i-1] = coefficients[i] * i

        return Polynomial(derivative_x)

    def y_derivative(self,polynomial:Polynomial) -> Polynomial:

        coefficients = polynomial.coefficients
        derivative_y = np.zeros(coefficients.shape)
        for i in range(1,coefficients.shape[1]):
            derivative_y[:,i-1] = coefficients[:,i] * i

        return Polynomial(derivative_y)

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
                          ) -> np.ndarray:
        # TODO: rewrite
        P = polynomials[0].coefficients
        Q = polynomials[1].coefficients

        # N = max(self.monomials[0].getPowers()[0], self.monomials[1].getPowers()[0])
        # M = max(self.monomials[0].getPowers()[1], self.monomials[1].getPowers()[1])

        # BIG_Matrix = np.zeros(shape=(4, P.shape[0] + N, P.shape[1] + M))

        T_monom1 = P * monomial_derivatives[0].coeff
        T_monom2 = Q * monomial_derivatives[1].coeff

        T_monom1 = Polynomial(T_monom1)
        T_monom2 = Polynomial(T_monom2)


        k1 = monomial_derivatives[0].powers[0]
        m1 = monomial_derivatives[0].powers[1]
        end11= k1 + P.shape[0]
        end12= m1 + P.shape[1]

        # BIG_Matrix[0, k1:end11, m1:end12] = T_monom1

        k2 = monomial_derivatives[1].powers[0]
        m2 = monomial_derivatives[1].powers[1]
        end21 = k2 + Q.shape[0]
        end22 = m2 + Q.shape[1]

        # BIG_Matrix[1, k2:end21, m2:end22] = T_monom2

        D_P_1 = self.monomials[0].coeff * polynomial_derivatives[0] * (-1)
        D_P_2 = self.monomials[1].coeff * polynomial_derivatives[1] * (-1)

        D_P_1 = Polynomial(D_P_1)
        D_P_2 = Polynomial(D_P_2)

        k3 = self.monomials[0].powers[0]
        m3 = self.monomials[0].powers[1]
        end31 = k3 + polynomial_derivatives[0].coefficients.shape[0]
        end32 = m3 + polynomial_derivatives[0].coefficients.shape[1]
        # BIG_Matrix[2, k3:end31, m3:end32] = D_P_1*(-1)

        k4 = self.monomials[1].powers[0]
        m4 = self.monomials[1].powers[1]
        end41 = k4 + polynomial_derivatives[1].coefficients.shape[0]
        end42 = m4 + polynomial_derivatives[1].coefficients.shape[1]
        # BIG_Matrix[3, k4:end41, m4:end42] = D_P_2*(-1)

        K1 = min(k1, k2, k3, k4)
        M1 = min(m1, m2, m3,m4)

        endK1 = max(end11, end21, end31, end41)
        endM1 = max(end12, end22, end32, end42)

        BIG_Matrix = [T_monom1,T_monom2,D_P_1,D_P_2]
        return BIG_Matrix,[K1,endK1,M1,endM1],[[k1,m1],[k2,m2],[k3,m3],[k4,m4]]

    def create_SOLE(self,polyPowers: list,
                    boundaries1: list,
                    boundaries2: list):
        cols = polyPowers[0]*polyPowers[1]*self.n
        rows = (boundaries1[1] - boundaries1[0] + 1)*(boundaries1[3] - boundaries1[2] + 1)
        rows += (boundaries2[1] - boundaries2[0] + 1)*(boundaries2[3] - boundaries2[2] + 1)

        return np.zeros((rows,cols))


    def fill_SOLE(self,BIG_Matrices: list,
                  shifts: list,
                  boundaries: list,
                  polyPowers: list,
                  SOLE: np.ndarray):
        # TODO: rewrite
        row_n = 0
        for m in range(len(BIG_Matrices)):
            BIG_MATRIX = BIG_Matrices[m]
            shift = shifts[m]
            boundary = boundaries[m]
            for i in range(boundary[0],boundary[1]):
                for j in range(boundary[2],boundary[3]):
                    row = np.zeros((1,SOLE.shape[1]))
                    for k in range(len(BIG_MATRIX)):
                        a_ij = BIG_MATRIX[k,i,j]
                        if k<2:
                            ind_i = i -shift[k][0]
                            ind_j = j -shift[k][1]
                        else:
                            ind_i = i - shift[k][0] + 1
                            ind_j = j - shift[k][1] + 1



                        position =  ind_i*polyPowers[1] + ind_j - 1
                        try:
                            row[0,position] += a_ij
                        except IndexError:
                            # print(k)
                            pass
                    SOLE[row_n,:] = row
                    row_n += 1



        return SOLE




    def commutator_search(self):
        [P,Q] = self.generateCommutator()
        print("P shape:",P.coefficients.shape)
        print("Q shape:",Q.coefficients.shape)
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
        print(shifts1)
        print(boundaries1)
        # print(BIG_Matrix1)

        print(shifts2)
        print(boundaries2)
        # print(BIG_Matrix2)
        SOLE = self.create_SOLE([P.shape[0],P.shape[1]],boundaries1,boundaries2)
        print("SOLE shape:",SOLE.shape)
        SOLE = self.fill_SOLE([BIG_Matrix1,BIG_Matrix2],[shifts1,shifts2],
                              [boundaries1,boundaries2],P.shape,SOLE)
        print(SOLE)

def __str__(self):

        return "{} d/dx + {} d/dy".format(self.monomials[0],self.monomials[1])



if __name__ == "__main__":
    component1 = Monomial(2,2,[2,3])
    component2 = Monomial(2,1,[3,4])
    moComm = MonomialCommutator(2,[component1,component2])
    moComm.commutator_search()



    # print(moComm)

    # polynomial = np.array([[1,1,2],[2,3,5]])
    # print("Polynomial:\n",polynomial)
    #
    # der_x = moComm.x_derivative(polynomial)
    # print("derivative x\n",der_x)
    #
    # der_y = moComm.y_derivative(polynomial)
    # print("derivative y\n",der_y)
    #
    # monomial = Monomial(3,2,[2,3,3])
    # print("Monomial:\n",monomial)
    # print("derivative x\n",moComm.monomial_derivative(monomial,1))
    # print("derivative y\n",moComm.monomial_derivative(monomial,2))
    # print("derivative z\n",moComm.monomial_derivative(monomial,3))
