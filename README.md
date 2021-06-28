# LBM-DEM-FEM-Coupling
The model enables heat transfer simulation of suspensions with ellipsoidal particles.

It's built based on following open-sorce codes:
-LBM: Palabos v2.1r0 (https://palabos.unige.ch/)
-DEM: LIGGGHTS v3.1.0 (https://www.cfdem.com/media/DEM/docu/Manual.html)
-LB-DEM Coupling: Philippe Seil (https://github.com/ParticulateFlow/LBDEMcoupling-public)
-FEM: QuickerSim on MATLAB® (https://quickersim.com/about/)

## Installation
1. install Palabos;
2. install LIGGGHTS:
    - git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git;
    - git checkout b3030a8 (to get the old version v3.1.0);
    - make yes-ASPHERE yes-RIGID (insert necessary packages, and copy necessary files from lammps/src);
    - copy files in the LIGGGHTS_Extra folder to LIGGGHTS-PUBLIC/src (these enable inserting pack of ellipsoids);
    - (make fedora; mpirun lmp_fedora < in.ellipsoid; run an example to see whether LIGGGHTS works well);
    - make makeshlib (update makefile.shlib file, what will be used to build the shared library);
    - make -f Makefile.shlib fedora_fpic (create shared library lammps.so);
3. install LBDEMCoupling:
    - git clone (this repository);
    - move the files "fix_lb_coupling_onetoone.cpp" and "fix_lb_coupling_onetoone.h" to your LIGGGHTS/src;
    - rebuild the LIGGGHTS library according to the step above;
4. install Matlab Runtime (same version as required by MatlabLibrary readme);
5. write soruce paths:
    - gedit ~/.bashrc (open bashrc file and add following paths to the file);
    - export PALABOS_ROOT=path/to/palabos/
    - export LIGGGHTS_ROOT=path/to/liggghts/
    - export LBDEM_ROOT=path/to/lbdem/
    - export LD_LIBRARY_PATH=path/to/liggghts/src/:path/to/the/simulating/case(because MatlabLibrary xxx.so and xxx.h is put in the folder):path/to/matlab/v98/runtime/glnxa64:path/to/matlab/v98/bin/glnxa64:path/to/matlab/v98/extern/bin/glnxa64:path/to/matlab/v98/sys/os/glnxa64:path/to/v98/sys/opengl/lib/glnxa64
    - source ~/.bashrc;
6. edit Makefile of each showcase:
    - projectFiles = (showcase).cpp ${LBDEM_ROOT}/src/liggghtsCouplingWrapper.cpp ${LBDEM_ROOT}/src/latticeDecomposition.cpp
    - libraryPaths = ${LIGGGHTS_ROOT}/src ${LBDEM_ROOT}/examples/ShearFlow /home/osboxes/matlab/v98/runtime/glnxa64 /home/osboxes/matlab/v98/bin/glnxa64
    - includePaths = ${LBDEM_ROOT}/src ${LIGGGHTS_ROOT}/src /home/osboxes/matlab/v98/extern/include
    - libraries    = liblammps.so libtestFunc.so libmwmclmcrrt.so libmwmclmcr.so

## Run the simulation
1. go to the showcase folder;
2. command "make" to compile the excutable file xxx.o;
3. excute the simulation with "./xxx".
