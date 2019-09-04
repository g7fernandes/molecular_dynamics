#!/bin/bash
read -p "[d]ebug with mpiifort, compile with mpi[i]fort or [c]ompile with OpenMPI+gfortran: " option

if [ $option == 'c' ] 
then
    echo "compiling with mpifort"
    /usr/lib64/mpi/gcc/openmpi/bin/mpifort mod1.f90 linkedlist.f90 mod0.f90 saida.f90 data.f90 matprint.f90 m_config.f90 randnormal.f90 lennard.f90 -o lennard.out
    echo "output: lennard.out"
    read -p "Run lennard.out? [y/n]" option2
    if [ $option2 = 'y' ]
    then
        /usr/lib64/mpi/gcc/openmpi/bin/mpirun -n 1 ./lennard.out
    fi
else
    if [ $option == "d" ]
    then
        echo "debuging with intel mpiifort"
        mpiifort mod1.f90 linkedlist.f90 mod0.f90 saida.f90 data.f90 matprint.f90 m_config.f90 randnormal.f90 lennard.f90 -debug -O0 -shared-intel -nocheck -fno-omit-frame-pointer -fasynchronous-unwind-tables -fexceptions -o lennard.out
        echo "output: lennard.out"
    else
        echo "compiling with intel mpiifort without optimization"
        mpiifort mod1.f90 linkedlist.f90 mod0.f90 saida.f90 data.f90 matprint.f90 m_config.f90 randnormal.f90 lennard.f90 -O0 -g -traceback -check all -o lennard.out
        echo "output: lennard.out"
    fi
    
fi




