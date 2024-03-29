# LBM-DEM-FEM-Coupling
The model enables heat transfer simulation of suspensions with ellipsoidal particles.

It's built based on following open-sorce codes:
-LBM: Palabos v2.1r0 (https://palabos.unige.ch/)
-DEM: LIGGGHTS v3.1.0 (https://www.cfdem.com/media/DEM/docu/Manual.html)
-LB-DEM Coupling: Philippe Seil (https://github.com/ParticulateFlow/LBDEMcoupling-public)
-FEM: QuickerSim on MATLAB® (https://quickersim.com/about/)

## Installation
1. install Palabos;
    - git clone https://gitlab.com/unigespc/palabos.git;
    - git checkout d1cfe102 (to get the old version v2.1r0, new version deleted scons compiler);
    - copy the codes in palabos_changed_code into src\complexDynamics folder;
2. install LIGGGHTS:
    - git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git;
    - git checkout b3030a8 (to get the old version v3.1.0);
    - make yes-ASPHERE yes-RIGID yes-SRD (insert necessary packages, and copy necessary files from lammps/src);
    - copy files in the LIGGGHTS_Extra folder to LIGGGHTS-PUBLIC/src (these enable inserting pack of ellipsoids);
    - (make fedora; mpirun lmp_fedora < in.ellipsoid; run an example to see whether LIGGGHTS works well);
    - make makeshlib (update makefile.shlib file, what will be used to build the shared library);
    - make -f Makefile.shlib fedora_fpic (create shared library lammps.so);
3. install LBDEMCoupling:
    - git clone (LB-DEM_Coupling repository);
    - move the files "fix_lb_coupling_onetoone.cpp" and "fix_lb_coupling_onetoone.h" to your LIGGGHTS/src;
    - rebuild the LIGGGHTS library according to the step above;
4. install Matlab Runtime (same version as required by MatlabLibrary readme);
5. write soruce paths:
    - gedit ~/.bashrc (open bashrc file and add following paths to the file);
    - export PALABOS_ROOT=path/to/palabos/
    - export LIGGGHTS_ROOT=path/to/liggghts/
    - export LBDEM_ROOT=path/to/lbdem/
    - export LD_LIBRARY_PATH=path/to/liggghts/src/:path/to/matlab/version/runtime/glnxa64:path/to/matlab/version/bin/glnxa64:path/to/matlab/version/extern/bin/glnxa64:path/to/matlab/version/sys/os/glnxa64:path/to/version/sys/opengl/lib/glnxa64
    - source ~/.bashrc;
6. edit Makefile of each showcase:
    - projectFiles = (showcase).cpp
                    ${LBDEM_ROOT}/src/liggghtsCouplingWrapper.cpp
                    ${LBDEM_ROOT}/src/latticeDecomposition.cpp
    - libraryPaths = ${LIGGGHTS_ROOT}/src
                    /path/to/showcase(for specific libMatlabFE.so)
                    /path/to/matlab/version/runtime/glnxa64
                    /path/to/matlab/version/bin/glnxa64
    - includePaths = ${LBDEM_ROOT}/src
                    ${LIGGGHTS_ROOT}/src
                    /path/to/matlab/version/extern/include
    - libraries    = liblammps.so libMatlabFE.so libmwmclmcrrt.so libmwmclmcr.so

## Run the simulation
1. go to the showcase folder;
2. command "make" to compile the object file xxx.o;
3. excute the simulation with "./xxx".


## References

Shu Q., Pan L., Kneer R., Rietz M., Rohlfs W., 2023: Effective thermal conductivity simulations of suspensions containing non-spherical particles in shear flow, International Journal of Heat and Mass Transfer, Band 204, 123808.

Shu Q., Rietz M., Kneer R., Rohlfs W., 2023: Thermal convection enhancement of laminar particle-laden flow in a square duct: fully resolved numerical investigation, Proceedings of the 17th International Heat Transfer Conference, Cape Town, Western Cape, South Africa.
