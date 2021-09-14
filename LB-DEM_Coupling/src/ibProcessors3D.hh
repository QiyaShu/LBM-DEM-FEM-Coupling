/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Original Author: Philippe Seil (philippe.seil@jku.at)
 * Adapted by Qiya Shu (shu@wsa.rwth-aachen.de)
 */

#include "ibDef.h"
#include "ibCompositeDynamics3D.h"
#include "ibDynamicsParticleData.h"
#include "utils.h"

#include "lammps.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

namespace plb{

  /*
   * implementation of SetSingleEllipsoid3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSingleEllipsoid3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    // modifies stuff, but we don't want updates running
    // after each particle set... should be modif::dataStructure
    // modified[0] = modif::dataStructure;
    modified[0] = modif::nothing; 
  }

  template<typename T, template<typename U> class Descriptor>
  SetSingleEllipsoid3D<T,Descriptor>* SetSingleEllipsoid3D<T,Descriptor>::clone() const 
  {
    return new SetSingleEllipsoid3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSingleEllipsoid3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          
          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;

          // this one actually helps, believe it or not
          __builtin_prefetch(particleData,1);
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

                    
          T const dx = xGlobal - x[0];
          T const dy = yGlobal - x[1];
          T const dz = zGlobal - x[2];

          T const dx_com = xGlobal - com[0];
          T const dy_com = yGlobal - com[1];
          T const dz_com = zGlobal - com[2];

          T const qth = acos(q[0])*2;
          T const qa = q[1]/sin(qth/2);
          T const qb = q[2]/sin(qth/2);
          T const qc = q[3]/sin(qth/2);


          T const sf = calcSolidFractionE(dx,dy,dz,s[0],s[1],s[2],qa,qb,qc,qth);

          T const sf_old = particleData->solidFraction;
          int const id_old = (int) particleData->partId;
          
          plint const decFlag = (sf > SOLFRAC_MIN) + 2*(sf_old > SOLFRAC_MIN);
          
          switch(decFlag){
          case 0: // sf == 0 && sf_old == 0
            setToZero(*particleData);
            break; // do nothing
          case 1: // sf > 0 && sf_old == 0
            setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break;
          case 2: // sf == 0 && sf_old > 0
            if( id_old == id ) // then particle has left this cell
              setToZero(*particleData);
            break; // else do nothing
          case 3: // sf > 0 && sf_old > 0
            if( sf > sf_old || id_old == id )
              setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break; // else do nothing
          }
          // if desired, initialize interior of ellipsoid with ellipsoid velocity
          if(initVelFlag && sf > SOLFRAC_MAX)
            cell.defineVelocity(particleData->uPart);

        }
      }
    }
    
  }

	
  template<typename T, template<typename U> class Descriptor>
  T SetSingleEllipsoid3D<T,Descriptor>::calcSolidFractionE(T const dx_, T const dy_, T const dz_, T const sa_, T const sb_, T const sc_, T const q0_, T const q1_, T const q2_, T const q3_)
  {
    // return calcSolidFractionRec(dx_,dy_,dz_,s_,q_);

    static plint const slicesPerDim = 5;
    static T const sliceWidth = 1./((T)slicesPerDim);
    static T const fraction = 1./((T)(slicesPerDim*slicesPerDim*slicesPerDim));
    
    // should be sqrt(3.)/2.
    // add a little to avoid roundoff errors
    static const T sqrt3half = (T) sqrt(3.1)/2.; 

    T const dx_o = dx_*(q0_*q0_*(1-cos(-q3_))+cos(-q3_))+dy_*(q0_*q1_*(1-cos(-q3_))-q2_*sin(-q3_))+dz_*(q0_*q2_*(1-cos(-q3_))+q1_*sin(-q3_));
    T const dy_o = dx_*(q1_*q0_*(1-cos(-q3_))+q2_*sin(-q3_))+dy_*(q1_*q1_*(1-cos(-q3_))+cos(-q3_))+dz_*(q1_*q2_*(1-cos(-q3_))-q0_*sin(-q3_));
    T const dz_o = dx_*(q2_*q0_*(1-cos(-q3_))-q1_*sin(-q3_))+dy_*(q2_*q1_*(1-cos(-q3_))+q0_*sin(-q3_))+dz_*(q2_*q2_*(1-cos(-q3_))+cos(-q3_));

    T const sa_p = sa_ + sqrt3half;
    T const sb_p = sb_ + sqrt3half;
    T const sc_p = sc_ + sqrt3half;

    T dist_p = (dx_o*dx_o)/(sa_p*sa_p) + (dy_o*dy_o)/(sb_p*sb_p) + (dz_o*dz_o)/(sc_p*sc_p);
    if (dist_p > 1) return 0;

    T const sa_m = sa_ - sqrt3half;
    T const sb_m = sb_ - sqrt3half;
    T const sc_m = sc_ - sqrt3half;

    T dist_m = (dx_o*dx_o)/(sa_m*sa_m) + (dy_o*dy_o)/(sb_m*sb_m) + (dz_o*dz_o)/(sc_m*sc_m);
    if (dist_m < 1) return 1;


    T dx_0[slicesPerDim],dy_0[slicesPerDim],dz_0[slicesPerDim];

    // pre-calculate d[xyz]_sq for efficiency
    for(plint i=0;i<slicesPerDim;i++){
      T const delta = -0.5 + ((T)i+0.5)*sliceWidth;
      T const dx = dx_+delta; 
      T const dy = dy_+delta; 
      T const dz = dz_+delta;
      dx_0[i] = dx*(q0_*q0_*(1-cos(-q3_))+cos(-q3_))+dy*(q0_*q1_*(1-cos(-q3_))-q2_*sin(-q3_))+dz*(q0_*q2_*(1-cos(-q3_))+q1_*sin(-q3_));
      dy_0[i] = dx*(q1_*q0_*(1-cos(-q3_))+q2_*sin(-q3_))+dy*(q1_*q1_*(1-cos(-q3_))+cos(-q3_))+dz*(q1_*q2_*(1-cos(-q3_))-q0_*sin(-q3_));
      dz_0[i] = dx*(q2_*q0_*(1-cos(-q3_))-q1_*sin(-q3_))+dy*(q2_*q1_*(1-cos(-q3_))+q0_*sin(-q3_))+dz*(q2_*q2_*(1-cos(-q3_))+cos(-q3_));
    }

    pluint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      for(plint j=0;j<slicesPerDim;j++){
        for(plint k=0;k<slicesPerDim;k++){
          T dist = (dx_0[i]*dx_0[i])/(sa_*sa_) + (dy_0[j]*dy_0[j])/(sb_*sb_) + (dz_0[k]*dz_0[k])/(sc_*sc_);
          n += (dist < 1);
        }
      }
    }

    return fraction*((T)n);
  }


  template<typename T, template<typename U> class Descriptor>
  void SetSingleEllipsoid3D<T,Descriptor>::setValues(IBdynamicsParticleData<T,Descriptor> &p,
                                                  T const sf, T const dx, T const dy, T const dz)
  {    
    p.uPart.from_cArray(v);
    if(omega != 0){
      p.uPart[0] += omega[1]*dz - omega[2]*dy;
      p.uPart[1] += -omega[0]*dz + omega[2]*dx; 
      p.uPart[2] += omega[0]*dy - omega[1]*dx;
    }
    p.solidFraction = sf;
    p.partId = id;
  }

	
  template<typename T, template<typename U> class Descriptor>
  void SetSingleEllipsoid3D<T,Descriptor>::setToZero(IBdynamicsParticleData<T,Descriptor> &p)
  {
    p.uPart[0] = 0;
    p.uPart[1] = 0;
    p.uPart[2] = 0;
    p.solidFraction = 0;
    p.partId = 0;
  }



  /*
   * implementation of SetSingleSphere3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    // modifies stuff, but we don't want updates running
    // after each particle set... should be modif::dataStructure
    // modified[0] = modif::dataStructure;
    modified[0] = modif::nothing; 
  }

  template<typename T, template<typename U> class Descriptor>
  SetSingleSphere3D<T,Descriptor>* SetSingleSphere3D<T,Descriptor>::clone() const 
  {
    return new SetSingleSphere3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          
          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;

          // this one actually helps, believe it or not
          __builtin_prefetch(particleData,1);
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

                    
          T const dx = xGlobal - x[0];
          T const dy = yGlobal - x[1];
          T const dz = zGlobal - x[2];

          T const dx_com = xGlobal - com[0];
          T const dy_com = yGlobal - com[1];
          T const dz_com = zGlobal - com[2];

          T const sf = calcSolidFraction(dx,dy,dz,r);

          T const sf_old = particleData->solidFraction;
          int const id_old = (int) particleData->partId;
          
          plint const decFlag = (sf > SOLFRAC_MIN) + 2*(sf_old > SOLFRAC_MIN);
          
          switch(decFlag){
          case 0: // sf == 0 && sf_old == 0
            setToZero(*particleData);
            break; // do nothing
          case 1: // sf > 0 && sf_old == 0
            setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break;
          case 2: // sf == 0 && sf_old > 0
            if( id_old == id ) // then particle has left this cell
              setToZero(*particleData);
            break; // else do nothing
          case 3: // sf > 0 && sf_old > 0
            if( sf > sf_old || id_old == id )
              setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break; // else do nothing
          }
          // if desired, initialize interior of sphere with sphere velocity
          if(initVelFlag && sf > SOLFRAC_MAX)
            cell.defineVelocity(particleData->uPart);

        }
      }
    }
    
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_)
  {
    // return calcSolidFractionRec(dx_,dy_,dz_,r_);

    static plint const slicesPerDim = 5;
    static T const sliceWidth = 1./((T)slicesPerDim);
    static T const fraction = 1./((T)(slicesPerDim*slicesPerDim*slicesPerDim));
    
    // should be sqrt(3.)/2.
    // add a little to avoid roundoff errors
    static const T sqrt3half = (T) sqrt(3.1)/2.; 

    T const dist = dx_*dx_ + dy_*dy_ + dz_*dz_;

    T const r_p = r_ + sqrt3half;
    if (dist > r_p*r_p) return 0;

    T const r_m = r_ - sqrt3half;
    if (dist < r_m*r_m) return 1;

    T const r_sq = r_*r_;
    T dx_sq[slicesPerDim],dy_sq[slicesPerDim],dz_sq[slicesPerDim];

    // pre-calculate d[xyz]_sq for efficiency
    for(plint i=0;i<slicesPerDim;i++){
      T const delta = -0.5 + ((T)i+0.5)*sliceWidth;
      T const dx = dx_+delta; dx_sq[i] = dx*dx;
      T const dy = dy_+delta; dy_sq[i] = dy*dy;
      T const dz = dz_+delta; dz_sq[i] = dz*dz;
    }

    pluint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      for(plint j=0;j<slicesPerDim;j++){
        for(plint k=0;k<slicesPerDim;k++){
          n += (dx_sq[i] + dy_sq[j] + dz_sq[k] < r_sq);
        }
      }
    }

    return fraction*((T)n);
  }



  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setValues(IBdynamicsParticleData<T,Descriptor> &p,
                                                  T const sf, T const dx, T const dy, T const dz)
  {    
    p.uPart.from_cArray(v);
    if(omega != 0){
      p.uPart[0] += omega[1]*dz - omega[2]*dy;
      p.uPart[1] += -omega[0]*dz + omega[2]*dx; 
      p.uPart[2] += omega[0]*dy - omega[1]*dx;
    }
    p.solidFraction = sf;
    p.partId = id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setToZero(IBdynamicsParticleData<T,Descriptor> &p)
  {
    p.uPart[0] = 0;
    p.uPart[1] = 0;
    p.uPart[2] = 0;
    p.solidFraction = 0;
    p.partId = 0;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x_,
                                                   T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_)
    : x(x_), force(force_), torque(torque_), wrapper(wrapper_)
  {}

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig)
    : BoxProcessingFunctional3D_L<T,Descriptor>(orig), x(orig.x),
      force(orig.force), torque(orig.torque), wrapper(orig.wrapper) {}
  
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();
    
    // "real" domain size is nx-2 etc
    plint nx = lattice.getNx()-2, ny = lattice.getNy()-2, nz = lattice.getNz()-2;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T,Descriptor>* particleData=
            getParticleDataFromCell<T,Descriptor>(cell);



          if(!particleData) continue;  //for the cell without particle, jump to the next

          // LIGGGHTS indices start at 1
          plint const id = particleData->partId;
          if(id < 1) continue; // no particle here

          plint const ind = wrapper.lmp->atom->map(id);
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

	  T dx = xGlobal - x[ind][0];  //x is the center position of the particle
	  T dy = yGlobal - x[ind][1];
	  T dz = zGlobal - x[ind][2];
          
          // minimum image convention, needed if
          // (1) Periodic BC are used and
          // (2) both ends of Periodic BC lie on the same processor
          if(dx > nx/2) dx -= nx;
          else if(dx < -nx/2) dx += nx;  //x direction periodic
          /*if(dy > ny/2) dy -= ny;
          else if(dy < -ny/2) dy += ny;
          if(dz > nz/2) dz -= nz;
          else if(dz < -nz/2) dz += nz;*/
	           
          T const forceX = particleData->hydrodynamicForce[0];
          T const forceY = particleData->hydrodynamicForce[1];
          T const forceZ = particleData->hydrodynamicForce[2];
          
          T const torqueX = dy*forceZ - dz*forceY;
          T const torqueY = -dx*forceZ + dz*forceX;
          T const torqueZ = dx*forceY - dy*forceX;

          addForce(ind,0,forceX);
          addForce(ind,1,forceY);
          addForce(ind,2,forceZ);
          
          addTorque(ind,0,torqueX);
          addTorque(ind,1,torqueY);
          addTorque(ind,2,torqueZ);
        }
      }
    }
  }
        

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint const partId, plint const coord, T const value)
  {
    force[3*partId+coord] += value;
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint const partId, plint const coord, T const value)
  {
    torque[3*partId+coord] += value;
  }

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>* SumForceTorque3D<T,Descriptor>::clone() const
  { 
    return new SumForceTorque3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
  }

  
};
