CentOS Linux
icc version 17.0.0 (gcc version 5.5.0 compatibility)

Compilation: icc -USE_OMP -qopenmp UJ_HH_Exp3_LC.cpp -o UJ_HH_Exp3_LC.out
Give permission to run: chmod +x run_script_g_syn_P_g_syn.sh
Run single sample: srun -N 1 -n 1 -p gpu 3000 ./run_script_g_syn_P_g_syn.sh