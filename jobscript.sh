#!/bin/bash

#mpirun --mca btl_openlib_warn_no_device_params_found 0 --mca btl ^openib -n 6 python3 matMul_V1_1.py
mpirun --mca btl_openlib_warn_no_device_params_found 0 --mca btl ^openib -n 6 python3 test1_mpi.py
#mpiexec --mca btl_openib_warn_no_device_params_found 0 --mca btl ^openib -n 6 python3 matMul_V2.py
