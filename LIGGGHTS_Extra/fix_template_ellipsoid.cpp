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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_ellipsoid.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsertEllipsoid.h"
#include "region.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "vector_liggghts.h"
#include "fix_region_variable.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;
using namespace FixConst;

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixTemplateEllipsoid::FixTemplateEllipsoid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (domain->dimension != 3)
    error->fix_error(FLERR,this,"this fix is for 3D simulations only");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->fix_error(FLERR,this,"not enough arguments");
  seed = atoi(arg[3]) + comm->me;
  random = new RanPark(lmp,seed);

  iarg = 4;

  // set default values
  atom_type = 2;
  vol_limit = 1e-14;

  pdf_shapex = NULL;
  pdf_shapey = NULL;
  pdf_shapez = NULL;
  pdf_density = NULL;

  delete pti;
  pti = new ParticleToInsertEllipsoid(lmp);

  n_pti_max = 0;
  pti_list = NULL;

  reg = NULL;
  reg_var = NULL;

  //parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"atom_type") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      atom_type=atoi(arg[iarg+1]);
      if (atom_type < 1) error->fix_error(FLERR,this,"invalid atom type (must be >=1)");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int ireg = domain->find_region(arg[iarg+1]);
      if (ireg < 0) error->fix_error(FLERR,this,"illegal region");
      reg = domain->regions[ireg];
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region_variable") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0) error->fix_error(FLERR,this,"illegal region/variable fix");
      reg_var = static_cast<FixRegionVariable*>(modify->fix[ifix]);
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"volume_limit") == 0)
    {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'volume_limit'");
      vol_limit = atof(arg[iarg+1]);
      if(vol_limit <= 0)
        error->fix_error(FLERR,this,"volume_limit > 0 required");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"shape") == 0)
    {
      pdf_shapex = new PDF(error);
      pdf_shapey = new PDF(error);
      pdf_shapez = new PDF(error);
      hasargs = true;
      if(strcmp(this->style,"particletemplate/ellipsoid"))
	error->fix_error(FLERR,this,"keyword shape only valid for particletemplate/ellipsoid");
      if (iarg+5 > narg)
        error->all(FLERR,"Illegal fix particletemplate/ellipsoid command, not enough arguments");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
        double values[] = { atof(arg[iarg+2])*force->cg(),//cg(atom_type),
                            atof(arg[iarg+3])*force->cg(),//cg(atom_type),
                            atof(arg[iarg+4])*force->cg()};//cg(atom_type);
        if( values[0] <= 0. or values[1] <= 0. or values[2] <= 0.)
          error->all(FLERR,"Illegal fix particletemplate/ellipsoid command, shape parameters must be >= 0");
        pdf_shapex->set_params<RANDOM_CONSTANT>(values[0]);
        pdf_shapey->set_params<RANDOM_CONSTANT>(values[1]);
        pdf_shapez->set_params<RANDOM_CONSTANT>(values[2]);
        iarg += 5;
      }
      else if (strcmp(arg[iarg+1],"gaussian") == 0)
      {
        if (iarg+6 > narg)
           error->fix_error(FLERR,this,"not enough arguments");
        double mux = atof(arg[iarg+2]);
        double muy = atof(arg[iarg+3]);
        double muz = atof(arg[iarg+4]);
        double sigma = atof(arg[iarg+5]);
        if(mux <= 0.)
          error->fix_error(FLERR,this,"illegal mux value for density");
        if(muy <= 0.)
          error->fix_error(FLERR,this,"illegal muy value for density");
        if(muz <= 0.)
          error->fix_error(FLERR,this,"illegal muz value for density");
        if( sigma <= 0. ) error->all(FLERR,"illegal sigmax value for density");
        pdf_shapex->set_params<RANDOM_GAUSSIAN>(mux,sigma);
        pdf_shapey->set_params<RANDOM_GAUSSIAN>(muy,sigma);
        pdf_shapez->set_params<RANDOM_GAUSSIAN>(muz,sigma);
        iarg += 6;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
        if (iarg+8 > narg)
           error->fix_error(FLERR,this,"not enough arguments");
        double shxmin = atof(arg[iarg+2]);
        double shxmax = atof(arg[iarg+3]);
        double shymin = atof(arg[iarg+4]);
        double shymax = atof(arg[iarg+5]);
        double shzmin = atof(arg[iarg+6]);
        double shzmax = atof(arg[iarg+7]);
        if(shxmax < shxmin)
          error->fix_error(FLERR,this,"max shapex less than min shapex");
        if(shymax < shymin)
          error->fix_error(FLERR,this,"max shapey less than min shapey");
        if(shzmax < shzmin)
          error->fix_error(FLERR,this,"max shapez less than min shapez");
        if(shxmin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapex");
        if(shymin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapey");
        if(shzmin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapez");

        pdf_shapex->set_params<RANDOM_UNIFORM>(shxmin,shxmax);
        pdf_shapey->set_params<RANDOM_UNIFORM>(shymin,shymax);
        pdf_shapez->set_params<RANDOM_UNIFORM>(shzmin,shzmax);
        iarg += 8;
      }
      else 
          error->fix_error(FLERR,this,"fix particletemplate/ellipsoid currently supports only constant, gaussian and uniform shape params");
     }
     else if (strcmp(arg[iarg],"shapex") == 0 or strcmp(arg[iarg],"shapey") == 0 or strcmp(arg[iarg],"shapez") == 0)
     {
        pdf_shapex = new PDF(error);
        pdf_shapey = new PDF(error);
        pdf_shapez = new PDF(error);
        hasargs = true;
        if(strcmp(this->style,"particletemplate/ellipsoid"))
          error->fix_error(FLERR,this,"keyword shape only valid for particletemplate/ellipsoid");
        if (iarg+3 > narg)
          error->all(FLERR,"Illegal fix particletemplate/ellipsoid command, not enough arguments");
        if (strcmp(arg[iarg+1],"constant") == 0)
        {
          double value = atof(arg[iarg+2])*force->cg();//cg(atom_type);
          if( value <= 0. )
            error->all(FLERR,"Illegal fix particletemplate/ellipsoid command, shape parameters must be >= 0");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_CONSTANT>(value);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_CONSTANT>(value);
          if(strcmp(arg[iarg],"shapez") == 0)
              pdf_shapez->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
        }
        else if (strcmp(arg[iarg+1],"gaussian") == 0)
        {
          if (iarg+4 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if(mu <= 0.)
            error->fix_error(FLERR,this,"illegal mux value for density");
          if( sigma <= 0. ) error->all(FLERR,"illegal sigma value for density");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          if(strcmp(arg[iarg],"shapez") == 0)
            pdf_shapez->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          iarg += 4;
        }
        else if (strcmp(arg[iarg+1],"uniform") == 0)
        {
          if (iarg+4 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          double shmin = atof(arg[iarg+2]);
          double shmax = atof(arg[iarg+3]);
          if(shmax < shmin)
            error->fix_error(FLERR,this,"max shape less than min shape");
          if(shmin <= 0.0)
            error->fix_error(FLERR,this,"Illegal value min shapex");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_UNIFORM>(shmin,shmax);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_UNIFORM>(shmin,shmax);
          if(strcmp(arg[iarg],"shapez") == 0)
            pdf_shapez->set_params<RANDOM_UNIFORM>(shmin,shmax);
          iarg += 4;
        }
        else
          error->fix_error(FLERR,this,"fix particletemplate/ellipsoid currently supports only constant, gaussian and uniform shape params");
    }
    else if (strcmp(arg[iarg],"density") == 0)
    {
      hasargs = true;
      if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments");
      pdf_density = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2]);
          if( value <= 0.) error->fix_error(FLERR,this,"density must be >= 0");
          pdf_density->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else
          error->fix_error(FLERR,this,"fix particletemplate/ellipsoid currently supports only constant density");
    }
     
    else if(strcmp(style,"particletemplate/ellipsoid") == 0)
        error->fix_error(FLERR,this,"unrecognized keyword");
  }

  if(pdf_density == NULL) error->fix_error(FLERR,this,"have to define 'density'");

  // end here for derived classes
  if(strcmp(this->style,"particletemplate/ellipsoid"))return;

  if(pdf_shapex == NULL) error->fix_error(FLERR,this,"have to define 'shapex'");
  if(pdf_shapey == NULL) error->fix_error(FLERR,this,"have to define 'shapey'");
  if(pdf_shapez == NULL) error->fix_error(FLERR,this,"have to define 'shapez'");

  // set mass and volume expectancy
  volume_expect = expectancy(pdf_shapex)*expectancy(pdf_shapey)*expectancy(pdf_shapez)*4.*M_PI/3.;
  mass_expect = expectancy(pdf_density) * volume_expect;


}

/* ---------------------------------------------------------------------- */

FixTemplateEllipsoid::~FixTemplateEllipsoid()
{
    delete random;

    delete pdf_density;
    delete pdf_shapex;
    delete pdf_shapey;
    delete pdf_shapez;

    if(strcmp(style,"particletemplate/ellipsoid") == 0)
    {
        delete pti;
        if(pti_list) delete_ptilist();
    }
}

/* ----------------------------------------------------------------------*/

int FixTemplateEllipsoid::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------*/

Region* FixTemplateEllipsoid::region()
{
    if(reg_var) return reg_var->region();
    else return reg;
}

/* ----------------------------------------------------------------------*/

void FixTemplateEllipsoid::randomize_single()
{
    
    pti->atom_type = atom_type;
    ParticleToInsertEllipsoid *pties_ptr = dynamic_cast<ParticleToInsertEllipsoid*>(pti);

    // randomize shape
    double shape[3];
    shape[0] = rand(pdf_shapex, random);
    shape[1] = rand(pdf_shapey, random);
    shape[2] = rand(pdf_shapez, random);
    pties_ptr->shape_ins[0] = shape[0];
    pties_ptr->shape_ins[1] = shape[1];
    pties_ptr->shape_ins[2] = shape[2];

    double radius_ = std::max(std::max(shape[0], shape[1]), shape[2]);
    //MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape, blockiness, &radius_);
    pti->radius_ins[0] = pti->r_bound_ins = radius_;

    // randomize density
    pti->density_ins = rand(pdf_density,random);

    // calculate volume and mass
    pti->volume_ins = shape[0] * shape[1] * shape[2] * 4.*M_PI/3.;
    pti->mass_ins = pti->density_ins*pti->volume_ins;

    // init insertion position
    vectorZeroize3D(pti->x_ins[0]);

    pti->groupbit = groupbit;

}

/* ----------------------------------------------------------------------*/

void FixTemplateEllipsoid::init_ptilist(int n_random_max)
{
    if(pti_list) error->one(FLERR,"invalid FixTemplateEllipsoid::init_list()");
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsertEllipsoid**) memory->smalloc(n_pti_max*sizeof(ParticleToInsertEllipsoid*),"pti_list");
    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsertEllipsoid(lmp);
}

/* ----------------------------------------------------------------------*/

void FixTemplateEllipsoid::delete_ptilist()
{
    if(n_pti_max == 0) return;

    for(int i = 0; i < n_pti_max; i++)
       delete pti_list[i];

    memory->sfree(pti_list);
    pti_list = NULL;
    n_pti_max = 0;
}

/* ----------------------------------------------------------------------*/

void FixTemplateEllipsoid::randomize_ptilist(int n_random,int distribution_groupbit)
{
    for(int i = 0; i < n_random; i++)
    {
        ParticleToInsertEllipsoid *pties_ptr = dynamic_cast<ParticleToInsertEllipsoid*>(pti_list[i]);
        pti_list[i]->atom_type = atom_type;

        // randomize shapes
        double shape[3];
        shape[0] = rand(pdf_shapex, random);
        shape[1] = rand(pdf_shapey, random);
        shape[2] = rand(pdf_shapez, random);
        pties_ptr->shape_ins[0] = shape[0];
        pties_ptr->shape_ins[1] = shape[1];
        pties_ptr->shape_ins[2] = shape[2];

        double radius_ = std::max(std::max(shape[0], shape[1]), shape[2]);
        //MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape, blockiness, &radius_);
        pti_list[i]->radius_ins[0] = pti_list[i]->r_bound_ins = radius_;


        // randomize density
        pti_list[i]->density_ins = rand(pdf_density,random);

        // calculate volume and mass
        pti_list[i]->volume_ins = shape[0] * shape[1] * shape[2] * 4.*M_PI/3.;
        pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;

        // init insertion position
        vectorZeroize3D(pti_list[i]->x_ins[0]);
        vectorZeroize3D(pti_list[i]->v_ins);
        vectorZeroize3D(pti_list[i]->omega_ins);

        pti_list[i]->groupbit = groupbit | distribution_groupbit; 

        //pti_list[i]->distorder = distorder;
    }
    
}

/* ----------------------------------------------------------------------*/

double FixTemplateEllipsoid::min_rad()
{
    
    return pdf_min(pdf_shapez);
}

/* ----------------------------------------------------------------------*/

double FixTemplateEllipsoid::max_rad()
{
    
    return pdf_max(pdf_shapex);
}

/* ----------------------------------------------------------------------*/

double FixTemplateEllipsoid::max_r_bound()
{
    return pdf_max(pdf_shapex);
}

/* ----------------------------------------------------------------------*/

double FixTemplateEllipsoid::volexpect()
{
    if(volume_expect < vol_limit)
    {
        
        error->fix_error(FLERR,this,"Volume expectancy too small. Change 'volume_limit' "
        "if you are sure you know what you're doing");
    }
    return volume_expect;
}

/* ----------------------------------------------------------------------*/

double FixTemplateEllipsoid::massexpect()
{
    return mass_expect;
}

/* ----------------------------------------------------------------------*/

int FixTemplateEllipsoid::number_ellipsoids()
{
    return 1;
}

/* ----------------------------------------------------------------------*/

int FixTemplateEllipsoid::maxtype()
{
    return atom_type;
}

/* ----------------------------------------------------------------------*/

int FixTemplateEllipsoid::mintype()
{
    return atom_type;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTemplateEllipsoid::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTemplateEllipsoid::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]) + comm->me;

  random->reset(seed);
}
