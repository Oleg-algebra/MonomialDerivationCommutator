#!/bin/bash

mpirun --mca btl_openlib_warn_no_device_params_found 0 --mca btl ^openib -n 6 python3 test_mpi.py --case 111 --it 100
#play -n synth 1 sin 440
