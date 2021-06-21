/* ----------------------------------------------------------------------
   Adapted from LIGGGHTS v3.8.0 particle type Superquadric to type Ellipsoid for v3.1.0.
   
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_PARTICLE_TO_INSERT_ELLIPSOID_H
#define LMP_PARTICLE_TO_INSERT_ELLIPSOID_H

#include "memory.h"
#include "particleToInsert.h"
#include "pointers.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsertEllipsoid : public ParticleToInsert
    {
     public:

        ParticleToInsertEllipsoid(LAMMPS* lmp,int ns = 1);

        virtual ~ParticleToInsertEllipsoid();


        // per-sphere shape, position
        // if atom_type_vector exists, each sphere has different type
        double shape_ins[3];
        double inertia_ins[3];
        double quat_ins[4];


        // velocity and omega at insertion
        
        //double v_ins[3];
        //double omega_ins[3];

        virtual int insert();
        virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        virtual int set_x_v_omega(double *,double *,double *,double *);

        virtual void scale_pti(double r_scale);
    };

}

#endif
