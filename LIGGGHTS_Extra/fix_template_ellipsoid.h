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

#ifdef FIX_CLASS

FixStyle(particletemplate/ellipsoid,FixTemplateEllipsoid)

#else

#ifndef LMP_FIX_TEMPLATE_ELLIPSOID_H
#define LMP_FIX_TEMPLATE_ELLIPSOID_H

#include "fix.h"
#include "probability_distribution.h"
#include "fix_template_sphere.h"
//#include "math_extra_liggghts_superquadric.h"


namespace LAMMPS_NS {

class FixTemplateEllipsoid : public Fix {
 public:

  FixTemplateEllipsoid(class LAMMPS *, int, char **);
  ~FixTemplateEllipsoid();

  // inherited from Fix
  virtual void post_create(){}
  virtual int setmask();
  void write_restart(FILE *);
  void restart(char *);

  // access to protected properties
  virtual double volexpect();           
  virtual double massexpect();          
  virtual double min_rad();
  virtual double max_rad();
  virtual double max_r_bound();
  virtual int number_ellipsoids();
  virtual int maxtype();
  virtual int mintype();
  class Region *region();

  // single particle generation, used by fix pour/dev
  virtual void randomize_single();    
  class ParticleToInsertEllipsoid *pti;

  // many particle generation, used by fix insert commands
  virtual void init_ptilist(int);
  virtual void delete_ptilist();
  virtual void randomize_ptilist(int,int);
  //virtual void direct_set_ptlist(const int i, const void * const data, const int distribution_groupbit, const int distorder);
  int n_pti_max;
  class ParticleToInsertEllipsoid **pti_list;

  virtual void finalize_insertion() {}

 protected:

  int iarg;

  class Region *reg;
  class FixRegionVariable *reg_var;

  // random generator
  class RanPark *random;
  int seed;

  // properties of particle template
  int atom_type;
  class LMP_PROBABILITY_NS::PDF *pdf_shapex;
  class LMP_PROBABILITY_NS::PDF *pdf_shapey;
  class LMP_PROBABILITY_NS::PDF *pdf_shapez;   
  class LMP_PROBABILITY_NS::PDF *pdf_density;
  double volume_expect;
  double mass_expect;
  double vol_limit;
};

}

#endif // LMP_FIX_TEMPLATE_ELLIPSOID_H
#endif // FIX_CLASS
