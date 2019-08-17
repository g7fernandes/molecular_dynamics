#!/bin/bash
mpiifort mod1.f90 linkedlist.f90 mod0.f90 saida.f90 data.f90 matprint.f90 m_config.f90 randnormal.f90 lennard.f90 -debug -O0 -shared-intel -nocheck -fno-omit-frame-pointer -fasynchronous-unwind-tables -fexceptions -o lennard.out
