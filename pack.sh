#!/bin/bash
read -e -p "Enter particle generator file or position files: " VAR1

if [ "$VAR1" = "" ]; then
    echo "No particle generator given."
    tar -cvf lennard_files.tar potential_analysis.py pos_e_vel_inicial_usando_sim_anterior.py csv2vtk_particles.py lennard settings.ini verify_settings.py visualizar.py
else
    tar -cvf lennard_files.tar pos_e_vel_inicial_usando_sim_anterior.py csv2vtk_particles.py lennard $VAR1 settings.ini verify_settings.py visualizar.py
fi

