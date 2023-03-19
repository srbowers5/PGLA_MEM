#!/bin/bash


#first_tr = int(sys.argv[1])
#num_tr = int(sys.argv[2])
#first_rep = int(sys.argv[3])
#num_rep = int(sys.argv[4])
#first_str = int(sys.argv[5])
#num_str = int(sys.argv[6])    # Must match interm program numstr
#size_run = int(sys.argv[7])
#intermPath = sys.argv[8]
#first_run = int(sys.argv[9])
#code_template = sys.argv[10]
#lib_dir = sys.argv[11]
#num_runs = num_str / size_run


first_tr=1
last_tr=1
first_rep=1
num_rep=12
first_step=150001
num_step=50000
step_per_run=5000
dir=DO_REST_PCPG1/
run_num=30
code_version=intermE2_pcpg_1.f90
extra_lib_path=pyth_utils

./do_interm.py $first_tr $last_tr $first_rep $num_rep $first_step $num_step $step_per_run $dir $run_num $code_version $extra_lib_path

