from mpi4py import MPI
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix
from random_matrix import rand_matrix
from time import time

comm = MPI.COMM_WORLD
size = comm.size
rank = comm.rank
matrix = None
result_vector = None
eps = 1e-3

start_global = None
local_time = 0
times = None
if rank == 0:
    times = np.empty(size, dtype=np.float64)
    start_global = time()
    print("MPI process started.")
    N = np.random.randint(1000, 1001)
    non_zeros = np.random.randint(N)
    # matrix = coo_matrix(rand_matrix(N, non_zeros))
    # matrix = mmread("m_t1/m_t1.mtx")
    # matrix = mmread("rail_79841/rail_79841.mtx")
    matrix = mmread("sparsine/sparsine.mtx")
    vector = coo_matrix(np.ones((1, matrix.shape[0])))

    matrix = matrix.tocsr()

    print(f'matrix shape: {matrix.shape}')
    print(f'values number: {matrix.nnz}')
    chunk = matrix.shape[0] // (size - 1)
    resid = matrix.shape[0] % (size - 1)
    print(f'chunk size: {chunk}')
    print(f"resid: {resid}")
    if resid != 0:
        from_ind1 = np.arange(0, resid * (chunk + 1), chunk + 1)
        to_ind1 = from_ind1 + (chunk + 1)
        from_ind2 = np.arange(to_ind1[-1], matrix.shape[0], chunk)
        to_ind2 = from_ind2 + chunk

        from_ind = np.concatenate((from_ind1, from_ind2))
        to_ind = np.concatenate((to_ind1, to_ind2))
    else:
        from_ind = np.arange(0, matrix.shape[0], chunk)
        to_ind = from_ind + chunk
    print(f"rank {rank} distributing data between other ranks")
    for i in range(size - 1):
        rows = []
        for j in range(from_ind[i], to_ind[i]):
            rows.append(matrix.getrow(j))
        comm.send([rows, vector], dest=i + 1)
        # print(f"rank {i+1} recived rows {to_ind[i] - from_ind[i]}")
    print("Distribution finished!!!")
    end_distrib = time()
    print(f'distribution_time: {end_distrib - start_global}')
    # comm.Barrier()
else:
    # comm.Barrier()
    # print(f'rank {rank} computing result....')

    rows, vect = comm.recv(source=0)
    start_local = time()
    # comm.Barrier()
    # print(f'rank {rank} computing result...')
    vect = coo_matrix(vect)
    result_vector = np.zeros(len(rows))
    for i in range(len(rows)):
        res = vect.multiply(rows[i]).sum()

        result_vector[i] = res
    # print(f'rank {rank} finished calculations')
    end_local = time()
    local_time = end_local - start_local
    comm.send([result_vector, local_time], dest=0)

    # print(f'rank {rank} --- computation_time: {local_time} ')

# comm.Barrier()

# comm.Gather(local_time,times,root=0)
if rank == 0:
    print('gathering data....')
    times = np.zeros(size - 1)
    # print(times)
    result_vector = np.asarray([])
    for i in range(1, size):
        coordinates, time_local = comm.recv(source=i)
        times[i - 1] = time_local
        # print(coordinates)
        result_vector = np.concatenate((result_vector, coordinates))
    print(f"max execution time: {times.max()}")
    print(f"min execution time: {times.min()}")
    print(f'average execution time: {times.sum() / (size - 1)}')
    isCorrect = True
    python_vect = matrix.dot(vector.transpose()).todense()
    for i in range(len(python_vect)):
        if (abs(python_vect[i] - result_vector[i]) > eps):
            print(f'position: {i}')
            print(f'python value: {python_vect[i]}')
            print(f'my value: {result_vector[i]}')
            isCorrect = False
            break
    message = "Status: "
    if not isCorrect:
        # print(result_vector == matrix.dot(vector.transpose()).todense())
        # print(f"my computations: {result_vector}")
        # print(f'python: {matrix.dot(vector.transpose()).todense()}')
        print(message + "Fail")
    else:
        print(message + "OK")

    print(f'total time: {time() - start_global}')


