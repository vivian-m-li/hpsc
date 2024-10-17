# 8/29

- on alpine, need to do module load mpi/gcc/etc
- mpicxx main.cpp -o main.exe --> outputs main.exe
- ls -lrt --> sorts files by reverse order
- mpirun -n 4 ./main.exe --> to run in parallel (on 4 processors)
  - mpirun -n 4 ./main.exe | grep 'myPE: 0' --> only look at process 0
- srun on alpine
