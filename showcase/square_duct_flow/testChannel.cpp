/*
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

/*
 * This example consists of a square channel that is inserted with 
 * non-spherical particles.  
 * Pressure difference boundary conditions are applied.
 * Heat flux by a certain temperature at the wall is computed.
 */

#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "periodicPressureFunctionals3D.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"
#include "atom_vec_ellipsoid.h"

#include "libMatlabFE.h"
#include "mclmcrrt.h"  
#include "mclmcr.h"  
#include "matrix.h"  
#include "mclcppclass.h"


using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor//descriptors::D3Q19Descriptor
#define BASEDYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())  //BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS IBcompositeDynamics<T, DESCRIPTOR>( new BASEDYNAMICS )
#define TEMPERATURE_DESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor


/*
 * Define output functions: fluid vtk output, particulate vtk output, matlab output
 */

template<class BlockLatticeT, class TemperatureBlockLatticeT>
void writeVTK(BlockLatticeT& lattice, TemperatureBlockLatticeT& Temperature_lattice, IncomprFlowParam<T> const& parameters, PhysUnits3D<T> const& units, plint iter){

	//VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), units.getPhysLength(1));
	ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 3, units.getPhysLength(1));

	vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));
	//vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", units.getPhysVel(1));
	//vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1.0/units.getPhysTime(1));

	std::unique_ptr<MultiScalarField3D<T>> temperature = computeDensity(Temperature_lattice);
	vtkOut.writeData<float>(*temperature, "temperature", 1.0);

	T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;
	MultiScalarField3D<T> p(*computeDensity(lattice)); // p=rho/3
	subtractInPlace(p,1.); //p=p-1 on every lattice inside the domain
	vtkOut.writeData<float>(p,"pressure",p_fact );

	MultiScalarField3D<T> tmp(lattice);
	IBscalarQuantity sf = SolidFraction;
	applyProcessingFunctional(new GetScalarQuantityFromDynamicsFunctional<T,DESCRIPTOR,T>(sf),
                            lattice.getBoundingBox(),lattice,tmp);
	vtkOut.writeData<float>(tmp,"solidfraction",1. ); 
}


template<class BlockLatticeT, class TemperatureBlockLatticeT>
void writechannelheatflux (BlockLatticeT& lattice, TemperatureBlockLatticeT& Temperature_lattice,
		      T D, T rho_f, T cp_f, T k_f, plint iter, PhysUnits3D<T> const& units,
		      const char *topfilename, const char *bottomfilename, const char *frontfilename, const char *backfilename, const char *convectionfilename){

	
	std::vector<T> topheatGradient;
	for (plint dx=0; dx <= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-2; dz++){
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample1 = Temperature_lattice->get(dx,Temperature_lattice->getBoundingBox().y1,dz);
			T T_topsample1 = cell_sample1.getDynamics().computeDensity(cell_sample1);
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample2 = Temperature_lattice->get(dx,Temperature_lattice->getBoundingBox().y1-1,dz);
			T T_topsample2 = cell_sample2.getDynamics().computeDensity(cell_sample2);
			T top_Gradient = (T_topsample1 - T_topsample2) / units.getPhysLength(1) ;
			topheatGradient.push_back(top_Gradient);
		}
	}
	T averagetopheatflux =(std::accumulate(std::begin(topheatGradient),std::end(topheatGradient),0.0))/topheatGradient.size() * k_f;

	ofstream topfluxoutfile;
	topfluxoutfile.open(topfilename,ios::app);
	topfluxoutfile << iter << "\t" << averagetopheatflux <<"\t"<< "\n";
	topfluxoutfile.close();

	
	std::vector<T> bottomheatGradient;
	for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dz =1; dz <= Temperature_lattice->getBoundingBox().z1-1; dz++){
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample1 = Temperature_lattice->get(dx,Temperature_lattice->getBoundingBox().y0+1,dz);
		 	T T_bottomsample1 = cell_sample1.getDynamics().computeDensity(cell_sample1);
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample2 = Temperature_lattice->get(dx,Temperature_lattice->getBoundingBox().y0,dz);
			T T_bottomsample2 = cell_sample2.getDynamics().computeDensity(cell_sample2);
			T bottom_Gradient = (T_bottomsample1 - T_bottomsample2) / units.getPhysLength(1) ;
			bottomheatGradient.push_back(bottom_Gradient);
		}
	}        
	T averagebottomheatflux = (std::accumulate(std::begin(bottomheatGradient),std::end(bottomheatGradient),0.0))/bottomheatGradient.size() * k_f;

	ofstream bottomfluxoutfile;
	bottomfluxoutfile.open(bottomfilename,ios::app);
	bottomfluxoutfile << iter << "\t" << averagebottomheatflux <<"\t"<< "\n";
	bottomfluxoutfile.close();


	std::vector<T> frontheatGradient;
	for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dy =1; dy <= Temperature_lattice->getBoundingBox().y1-1; dy++){
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample1 = Temperature_lattice->get(dx, dy, Temperature_lattice->getBoundingBox().z1);
		 	T T_frontsample1 = cell_sample1.getDynamics().computeDensity(cell_sample1);
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample2 = Temperature_lattice->get(dx, dy, Temperature_lattice->getBoundingBox().z1-1);
			T T_frontsample2 = cell_sample2.getDynamics().computeDensity(cell_sample2);
			T front_Gradient = (T_frontsample1 - T_frontsample2) / units.getPhysLength(1) ;
			frontheatGradient.push_back(front_Gradient);
		}
	}        
	T averagefrontheatflux = (std::accumulate(std::begin(frontheatGradient),std::end(frontheatGradient),0.0))/frontheatGradient.size() * k_f;

	ofstream frontfluxoutfile;
	frontfluxoutfile.open(frontfilename,ios::app);
	frontfluxoutfile << iter << "\t" << averagefrontheatflux <<"\t"<< "\n";
	frontfluxoutfile.close();


	std::vector<T> backheatGradient;
	for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dy =0; dy <= Temperature_lattice->getBoundingBox().y1-2; dy++){
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample1 = Temperature_lattice->get(dx, dy, Temperature_lattice->getBoundingBox().z0+1);
		 	T T_backsample1 = cell_sample1.getDynamics().computeDensity(cell_sample1);
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample2 = Temperature_lattice->get(dx, dy, Temperature_lattice->getBoundingBox().z0);
			T T_backsample2 = cell_sample2.getDynamics().computeDensity(cell_sample2);
			T back_Gradient = (T_backsample1 - T_backsample2) / units.getPhysLength(1) ;
			backheatGradient.push_back(back_Gradient);
		}
	}        
	T averagebackheatflux = (std::accumulate(std::begin(backheatGradient),std::end(backheatGradient),0.0))/backheatGradient.size() * k_f;

	ofstream backfluxoutfile;
	backfluxoutfile.open(backfilename,ios::app);
	backfluxoutfile << iter << "\t" << averagebackheatflux <<"\t"<< "\n";
	backfluxoutfile.close();


	T averageheatflux=(abs(averagetopheatflux)+abs(averagebottomheatflux)+abs(averagefrontheatflux)+abs(averagebackheatflux))/4;
		
	
	/*
	T sumTemperature=0;
	for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dy =0; dy <= Temperature_lattice->getBoundingBox().y1-1; dy++){
			for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-1; dz++){
				Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_point = Temperature_lattice->get(dx, dy, dz);
				sumTemperature += cell_point.getDynamics().computeDensity(cell_point);
			}
		}
	}
	T averageTemperature=sumTemperature/(Temperature_lattice->getBoundingBox().x1*Temperature_lattice->getBoundingBox().y1*Temperature_lattice->getBoundingBox().z1);
	*/
	T sumBulkTemperature1=0;
	T sumVelocity1=0;
	for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dy =0; dy <= Temperature_lattice->getBoundingBox().y1-1; dy++){
			for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-1; dz++){
				Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_temp = Temperature_lattice->get(dx, dy, dz);
				Cell<T,DESCRIPTOR> const& cell_vel = lattice->get(dx, dy, dz);
				Array<T, DESCRIPTOR<T>::d> u1;
				cell_vel.getDynamics().computeVelocity(cell_vel,u1);
				sumBulkTemperature1 += cell_temp.getDynamics().computeDensity(cell_temp)*u1[0];
				sumVelocity1 += u1[0];
			}
		}
	}
	T bulkTemperature1=sumBulkTemperature1/sumVelocity1;
	
	/*
	T sumBulkTemperature2=0;
	T sumVelocity2=0;
	//for (plint dx=0; dx<= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dy =0; dy <= Temperature_lattice->getBoundingBox().y1-1; dy++){
			for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-1; dz++){
				Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_temp = Temperature_lattice->get(Temperature_lattice->getBoundingBox().x1-2, dy, dz);
				Cell<T,DESCRIPTOR> const& cell_vel = lattice->get(Temperature_lattice->getBoundingBox().x1-2, dy, dz);
				Array<T, DESCRIPTOR<T>::d> u2;
				cell_vel.getDynamics().computeVelocity(cell_vel,u2);
				sumBulkTemperature2 += cell_temp.getDynamics().computeDensity(cell_temp)*u2[0];
				sumVelocity2 += u2[0];
			}
		}
	//}
	T bulkTemperature2=sumBulkTemperature2/sumVelocity2;
	*/
	
	T averageConvectionCoefficient=averageheatflux/(1-bulkTemperature1);
	T averageNu=averageConvectionCoefficient*D/k_f;
	//T averageNu=rho_f*cp_f*(sumVelocity1/Temperature_lattice->getBoundingBox().y1/Temperature_lattice->getBoundingBox().z1)*D*D*(bulkTemperature1-bulkTemperature2)/units.getPhysLength(1)/(1-bulkTemperature1)/k_f/4;

	ofstream convectionoutfile;
	convectionoutfile.open(convectionfilename,ios::app);
	convectionoutfile << iter << "\t" << averageheatflux <<"\t"<< bulkTemperature1 <<"\t"<< averageConvectionCoefficient <<"\t"<<averageNu <<"\t"<< "\n";
	//convectionoutfile << iter << "\t" << sumVelocity2/Temperature_lattice->getBoundingBox().y1/Temperature_lattice->getBoundingBox().z1 << "\t" << bulkTemperature2 << "\t" << sumVelocity1/Temperature_lattice->getBoundingBox().y1/Temperature_lattice->getBoundingBox().z1 <<"\t"<< bulkTemperature1 <<"\t"<<averageNu <<"\t"<< "\n";
	convectionoutfile.close();


	pcout<<"writing the heat flux and convection coefficient data" << std::endl;

 }




template<class TemperatureBlockLatticeT>
T TrilinearInterpolation(TemperatureBlockLatticeT& TempLattice, Array<T,3> OuterVertices){
	T OuterVerticesTemperature;
	
	T dx = OuterVertices[0]-(int)OuterVertices[0];
	T dy = OuterVertices[1]-(int)OuterVertices[1];
	T dz = OuterVertices[2]-(int)OuterVertices[2];

//pcout << "OuterVertices: " << OuterVertices[0]<< ", " << OuterVertices[1] << ", " << OuterVertices[2]  << std::endl; 


        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell000 = TempLattice->get((int)OuterVertices[0],(int)OuterVertices[1],(int)OuterVertices[2]);
        T T000 = cell000.getDynamics().computeDensity(cell000);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell001 = TempLattice->get((int)OuterVertices[0],(int)OuterVertices[1],(int)OuterVertices[2]+1);
        T T001 = cell001.getDynamics().computeDensity(cell001);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell010 = TempLattice->get((int)OuterVertices[0],(int)OuterVertices[1]+1,(int)OuterVertices[2]);
        T T010 = cell010.getDynamics().computeDensity(cell010);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell011 = TempLattice->get((int)OuterVertices[0],(int)OuterVertices[1]+1,(int)OuterVertices[2]+1);
        T T011 = cell011.getDynamics().computeDensity(cell011);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell100 = TempLattice->get((int)OuterVertices[0]+1,(int)OuterVertices[1],(int)OuterVertices[2]);
        T T100 = cell100.getDynamics().computeDensity(cell100);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell101 = TempLattice->get((int)OuterVertices[0]+1,(int)OuterVertices[1],(int)OuterVertices[2]+1);
        T T101 = cell101.getDynamics().computeDensity(cell101);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell110 = TempLattice->get((int)OuterVertices[0]+1,(int)OuterVertices[1]+1,(int)OuterVertices[2]);
        T T110 = cell110.getDynamics().computeDensity(cell110);
        Cell<T,TEMPERATURE_DESCRIPTOR> const& cell111 = TempLattice->get((int)OuterVertices[0]+1,(int)OuterVertices[1]+1,(int)OuterVertices[2]+1);
        T T111 = cell111.getDynamics().computeDensity(cell111);

	T c0 = T000;
	T c1 = T100-T000;
	T c2 = T010-T000;
	T c3 = T001-T000;
	T c4 = T110-T010-T100+T000;
	T c5 = T011-T001-T010+T000;
	T c6 = T101-T001-T100+T000;
	T c7 = T111-T011-T101-T110+T100+T001+T010-T000;

	OuterVerticesTemperature = c0+c1*dx+c2*dy+c3*dz+c4*dx*dy+c5*dy*dz+c6*dz*dx+c7*dx*dy*dz;


	return OuterVerticesTemperature;
}


class SurfaceTemperatureFunc {
public:
    SurfaceTemperatureFunc (mwArray t_)
        : t(t_)
    { }
    T operator()(plint const& globalVertexIds)
    {

        T temperature;
	temperature = t.Get(1,globalVertexIds+1);
	return temperature;
    }
private:
    mwArray t;
};



/*template <typename T>
class PoiseuilleDensityAndVelocity {
public:
    PoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, T uMax_)
        : parameters(parameters_),
          uMax(uMax_),
	  R((parameters.getNy()-1)/2)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = uMax*(1-(((T)iY-R)*((T)iY-R)+((T)iZ-R)*((T)iZ-R))/R/R);
        u[1] = 0.;
        u[2] = 0.;
    }
private:
    IncomprFlowParam<T> parameters;
    T uMax;
    T R;
};

template <typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, T uMax_)
        : parameters(parameters_),
          uMax(uMax_),
	  R((parameters.getNy()-1)/2)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = uMax*(1-(((T)iY-R)*((T)iY-R)+((T)iZ-R)*((T)iZ-R))/R/R);
        u[1] = 0.;
        u[2] = 0.;
    }
private:
    IncomprFlowParam<T> parameters;
    T uMax;
    T R;
};*/




/*
 * Define parameters and transform into LB units
 */
int main(int argc, char* argv[]) {

	plbInit(&argc, &argv);

	// setting output directories
	std::string lbOutDir("results/"), demOutDir("results/");
	lbOutDir.append("tmp/"); demOutDir.append("post/");
	global::directories().setOutputDir(lbOutDir);


	char **argv_lmp = 0;
	argv_lmp = new char*[1];
	argv_lmp[0] = argv[0];

	LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

	const T rho_f = 997.1; //kg/m³
	const T nu_f = 1e-6; //m²/s
	const T cp_f = 4179; //J/(kgK)
	const T k_f = 0.613; //W/(mK)

	const double rho_s = 3970; //kg/m³
	const double cp_s = 765; //J/(kgK)
	const double k_s = 25; //W/(mK)

	T L = 4e-4; //duct width
	const T lx = 2*L, ly = L, lz = L;
	T d0 = 1e-4; //particle diameter
	T N = 100; //resolution

	T dp = 6.44; //pressure difference Pa
	T u0 = (dp/lx)*L*L*0.421731044868524/(12*nu_f*rho_f); //average velocity
	//T u_max = 2*u0;

	T dx = 1.0/N;
	T dt = 1.0/(N*N); // time resolution
	T u0_LB = dt/dx; 
	
	T t_0 = lx/u0; // reference time
	T dt_phys = t_0 * dt;
        T t_total = 1; // total simulated physical time space

	// this is equivalent to the variable command in LIGGGHTS/LAMMPS
	double r[3] = {d0/2, d0/4, d0/2};
	double k[3] = {k_s, k_s, k_s};//thermal conduction in different directions

	PhysUnits3D<T> units(L,u0,nu_f,lx,ly,lz,N,u0_LB,rho_f);

	IncomprFlowParam<T> parameters(units.getLbParam());


	T dpdx_LB = dp/lx/rho_f*(pow(L/u0,2)/L)*(pow(parameters.getDeltaT(),2)/parameters.getDeltaX());
pcout << "dp/dx_LB: " << dpdx_LB << std::endl;

	const plint maxSteps = ceil(t_total/dt_phys);
	const plint vtkSteps = 2000;
	//const plint gifSteps = 10;
	//const plint logSteps = 10;
 	int ibIter = 5; // Iterations for the immersed boundary method.	

	//temperatur parameters setting 
	T Prandtl = nu_f/(k_f/(rho_f*cp_f));
	//T LatticeDiff = parameters.getLatticeNu() /Prandtl;
	T LatticeDiff = (nu_f*units.getLbLength(1)*units.getLbLength(1)/units.getLbTime(1)) /Prandtl;//two options are the same

	T omega_temperature = (T) 1 / (TEMPERATURE_DESCRIPTOR<T>::invCs2 * LatticeDiff + (T) 0.5); //Relaxation parameter for the temperature
	
	pcout << " ---------------------------------------------- \n"
	      << "omega: " << parameters.getOmega() << "\n" //attention: omega between 0.6 and 0.9 for immersed boundary hydrodynamicForce (Noble1998)
	      << "omega_temperature: " << omega_temperature << "\n"
	      << "dt_phys: " << dt_phys << "\n"
	      << "average velocity: " << u0 << "\n"
	      << "Re : " << parameters.getRe() << "\n"
	      << "Pr : " << Prandtl << "\n"
	      << "vtkSteps: " << vtkSteps << "\n"
	      << "maxSteps: " << maxSteps << "\n"
	      << "grid size: " << parameters.getNx() << " " << parameters.getNy() << " " << parameters.getNz() << " \n"
	      << " ---------------------------------------------- " << std::endl;



/*
 * Setting up particles
 */

	// executes a LIGGGHTS input script
	wrapper.setVariable("xLength",lx);
	wrapper.setVariable("width",L);
	wrapper.setVariable("rx_part",r[0]);
	wrapper.setVariable("ry_part",r[1]);
	wrapper.setVariable("rz_part",r[2]);
	wrapper.setVariable("mass",rho_s*(4*M_PI*r[0]*r[1]*r[2]/3));
	wrapper.setVariable("density",rho_s);
	wrapper.setVariable("v_frac",0.05);//volume fraction of inserted particles

	wrapper.execFile("in.lbdem");


/*
 * Set up velocity lattice and temperature lattice for fluid
 */

	T omega = parameters.getOmega();
	/*T deltaRho = units.getLbRho(dp); //pressure difference in lattice units
	T rhoHi = 1., rhoLo = 1.-deltaRho;
pcout<< "deltaRho:" << deltaRho << std::endl;*/


	// get lattice decomposition from LIGGGHTS and create lattice according to parallelization
	// given in the LIGGGHTS input script
	LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),wrapper.lmp);
	SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
	ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
	plint envelopeWidth = 1;
 	plint largeEnvelopeWidth = 4; //Because of immersed walls. This defines the vicinty over which you will want to communicate things on your lattice. There are "halo" nodes that are used to communicate data in parallel executions


	pcout << "Generating multi-blocks." << std::endl;
	MultiBlockLattice3D<T, DESCRIPTOR> 
		*lattice = new MultiBlockLattice3D<T,DESCRIPTOR> (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
		defaultMultiBlockPolicy3D().getBlockCommunicator(),
		defaultMultiBlockPolicy3D().getCombinedStatistics(),
		defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
		new DYNAMICS );

	defineDynamics(*lattice,lattice->getBoundingBox(),new DYNAMICS);
	

	lattice->toggleInternalStatistics(false);

	/*
	MultiScalarField3D<T> *rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, largeEnvelopeWidth).release();
	rhoBar->toggleInternalStatistics(false);

	MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, largeEnvelopeWidth).release();
	j->toggleInternalStatistics(false);

	std::vector<MultiBlock3D*> rhoBarJarg;
	rhoBarJarg.push_back(lattice);
	rhoBarJarg.push_back(rhoBar);
	rhoBarJarg.push_back(j);
	integrateProcessingFunctional(
		new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(),
		lattice->getBoundingBox(), rhoBarJarg, 0);
	integrateProcessingFunctional(
		new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
		lattice->getBoundingBox(), rhoBarJarg, 2);*/


   	// temperature lattice 	   
	MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR > *Temperature_lattice = new MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR>(
		parameters.getNx(),parameters.getNy(),parameters.getNz(), 
		new AdvectionDiffusionRLBdynamics<T, TEMPERATURE_DESCRIPTOR>(omega_temperature) );
	Temperature_lattice->toggleInternalStatistics(false);
    
	MultiScalarField3D<T> *temperatureBar = generateMultiScalarField<T>((MultiBlock3D&) *Temperature_lattice, largeEnvelopeWidth).release();
	temperatureBar->toggleInternalStatistics(false);
   
	MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *Temperature_lattice, largeEnvelopeWidth).release();
	j->toggleInternalStatistics(false);



	//   Updating the body-force(gravity) in the whole domain in x-direction (equals to applying pressure difference)
        Array<T,DESCRIPTOR<T>::d> force(dpdx_LB,0.,0.);
	setExternalVector(*lattice,lattice->getBoundingBox(),
                          DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);

	lattice->periodicity().toggle(0,true);
	lattice->initialize();


	// Temperature boundary conditions
	const plint nx = parameters.getNx();
	const plint ny = parameters.getNy();
	const plint nz = parameters.getNz();
	Box3D topWall = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
	Box3D bottomWall = Box3D(0, nx-1, 0, 0, 0, nz-1);
	Box3D frontWall = Box3D(0, nx-1, 0, ny-1, nz-1, nz-1);
	Box3D backWall = Box3D(0, nx-1, 0, ny-1, 0, 0);
	Box3D everythingButWall = Box3D(0, nx-1, 1, ny-2, 1, nz-2);

	plint temperatureTimeFactor = 1;// Ratio between the temperature and the fluid time steps
	T T_wall = 1;//constant wall temperature
	T T_everythingButLid = 0;


	OnLatticeAdvectionDiffusionBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR> *bc = createLocalAdvectionDiffusionBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR>();
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, topWall, boundary::dirichlet);
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, bottomWall, boundary::dirichlet);
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, frontWall, boundary::dirichlet);
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, backWall, boundary::dirichlet);
    
	setBoundaryDensity(*Temperature_lattice, bottomWall, T_wall);
	setBoundaryDensity(*Temperature_lattice, topWall, T_wall);
	setBoundaryDensity(*Temperature_lattice, frontWall, T_wall);
	setBoundaryDensity(*Temperature_lattice, backWall, T_wall);
	Temperature_lattice->periodicity().toggle(0, true);
	Temperature_lattice->periodicity().toggle(1, false);
	Temperature_lattice->periodicity().toggle(2, false);

        temperatureBar->periodicity().toggle(0, true);
        temperatureBar->periodicity().toggle(1, false);
        temperatureBar->periodicity().toggle(2, false);

        j->periodicity().toggle(0, true);
        j->periodicity().toggle(1, false);
        j->periodicity().toggle(2, false);

	initializeAtEquilibrium(*Temperature_lattice, everythingButWall, T_everythingButLid, Array<T,3>((T)0,(T)0,(T)0));
	//initializeAtEquilibrium(*Temperature_lattice, topWall, T_wall, Array<T,3>((T)0,(T)0,(T)0));
	//initializeAtEquilibrium(*Temperature_lattice, bottomWall, T_wall, Array<T,3>((T)0,(T)0,(T)0));
	Temperature_lattice->initialize();

	integrateProcessingFunctional(new LatticeToPassiveAdvDiff3D<T,DESCRIPTOR,TEMPERATURE_DESCRIPTOR>((T)temperatureTimeFactor), lattice->getBoundingBox(),*lattice,*Temperature_lattice,1); //coupling between velocity lattice and temperature lattice
    
	bool incompressibleModel = true;

	std::vector<MultiBlock3D*> rhoBarJarg_temperature;
	rhoBarJarg_temperature.push_back(Temperature_lattice);
	rhoBarJarg_temperature.push_back(temperatureBar);
	rhoBarJarg_temperature.push_back(j);
	integrateProcessingFunctional(new ExternalRhoJcollideAndStream3D<T,TEMPERATURE_DESCRIPTOR>(), Temperature_lattice->getBoundingBox(), rhoBarJarg_temperature, 0);//RhoJ is copied back to lattice after collide and stream
	integrateProcessingFunctional(new BoxRhoBarJfunctional3D<T,TEMPERATURE_DESCRIPTOR>(), Temperature_lattice->getBoundingBox(), rhoBarJarg_temperature, 2);//lattice copied to RhoJ

	applyProcessingFunctional(new BoxRhoBarJfunctional3D<T,TEMPERATURE_DESCRIPTOR>(), Temperature_lattice-> getBoundingBox(), rhoBarJarg_temperature);

	// The next container block is necessary for the immersed-wall algorithm.
	MultiContainerBlock3D container(*temperatureBar);





/*	// Create the cylinder surface as a set of triangles.
	Array<T,3> originalCenter(0.0, parameters.getNy()/2, parameters.getNz()/2);
	TriangleSet<T> triangleSet;
	triangleSet = constructCylinder<T>(originalCenter, parameters.getNy(), parameters.getNz(), parameters.getNx(), 320, 100);

	// The next few lines of code are typical. They transform the surface geometry of the
	//   tube to more efficient data structures that are internally used by palabos.
	//   The TriangleBoundary3D structure will be later used to assign proper boundary conditions.
	DEFscaledMesh<T> defMesh(triangleSet, parameters.getNy(), yDirection, margin, extraLayer);
	defMesh.getMesh().inflate();
	TriangleBoundary3D<T> boundary(defMesh);


    // Voxelize the domain means: decide which lattice nodes are inside the solid
    // and which are outside.
    pcout << std::endl << "Voxelizing the simulation domain." << std::endl;
    int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain(triangleBoundary, flowType, param.bbox(), borderWidth,
            extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    {
        VtkImageOutput3D<T> vtkOut(outDir + "voxels_full_domain", param.dx);
        vtkOut.writeData<float>(*copyConvert<int,T>(voxelizedDomain.getVoxelMatrix(), param.bbox()), "voxel", 1.0);
    }*/


/*
 * Initializations
 */

	plint demSubsteps = 10;
	T dt_dem = dt_phys/(T)demSubsteps;

	// more variables to set: DEM timestep, number of steps to write a dumpfile
	// and dump directory
	wrapper.setVariable("t_step",dt_dem);
	wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
	wrapper.setVariable("dmp_dir",demOutDir);
	// executing second input file
	wrapper.execFile("in2.lbdem");
	// runs LIGGGHTS up to a certain step. Equivalent to the "run upto" command in LIGGGHTS/LAMMPS.
	wrapper.runUpto(demSubsteps-1);

	clock_t start = clock();
	clock_t loop = clock();
	clock_t end = clock(); 

	// Set up the application state for the MATLAB Runtime instance created in the application.
	if (!mclInitializeApplication(NULL,0)) {
		std::cerr << "could not initialize the application properly"<< std::endl;
		return -1;
	}
    
	// initialize lib  
	if( !libMatlabFEInitialize())  {  
		std::cout << "Could not initialize libmyFunc!" << std::endl;  
		return -1;  
	}

	//allocate space for variables for matlab
	plint nparticles = wrapper.lmp->atom->nlocal;
	mwArray SurfaceTemperatures(1530, 1, mxDOUBLE_CLASS);
	mwArray ParticleTemperatures(4357, 1, mxDOUBLE_CLASS);
	std::vector<T> ParticleTemperature(nparticles*ParticleTemperatures.NumberOfElements(),0); //set initial temperature for the paricle
	mwArray heatflux(1530, 1, mxDOUBLE_CLASS);// perpendicular heat flux of each node on the surface
	mwArray IT(1, 1, mxINT32_CLASS);
	mwArray VtkSteps(1, 1, mxINT32_CLASS);
	mwArray Dt_phys(1, 1, mxDOUBLE_CLASS);

	VtkSteps(1,1) = vtkSteps;
	Dt_phys(1,1) = dt_phys;


	//particleThermalConductionSetup(Radius, rho_s, cp_s, k_s);

	//TriangleSet<T> SphereMesh = TriangleSet<T>("sphereSurface.stl", (T) 1e-6);


	ofstream topfluxinfile;
	topfluxinfile.open("averagetopflux.txt",ios::trunc);
	topfluxinfile.close();
	ofstream bottomfluxoutfile;
	bottomfluxoutfile.open("averagebottomflux.txt",ios::trunc);
	bottomfluxoutfile.close();
	ofstream frontfluxinfile;
	frontfluxinfile.open("averagefrontflux.txt",ios::trunc);
	frontfluxinfile.close();
	ofstream backfluxoutfile;
	backfluxoutfile.open("averagebackflux.txt",ios::trunc);
	backfluxoutfile.close();


/*
 * Main loops
 */
	for (plint iT=0; iT<=maxSteps; ++iT) {//maxSteps

		static bool initWithVel = true;
		//setSpheresOnLattice(*lattice,wrapper,units,initWithVel);
		setEllipsoidOnLattice(*lattice,wrapper,units,initWithVel);
		if(initWithVel) initWithVel = false;

		if(iT%vtkSteps == 0 && iT > 0) { // LIGGGHTS does not write at timestep 0
			writeVTK(*lattice, *Temperature_lattice, parameters, units, iT);
			writechannelheatflux (lattice, Temperature_lattice, L, rho_f, cp_f, k_f, iT, units, "averagetopflux.txt", "averagebottomflux.txt", "averagefrontflux.txt", "averagebackflux.txt", "convectioncoefficient.txt");
		}




		/* Velocity lattice collision and stream, run liggghts */
		lattice->collideAndStream();


		getForcesFromLattice(*lattice,wrapper,units);
  
		// equilvalent to the "run" command in LIGGGHTS/LAMMPS
		wrapper.run(demSubsteps);

		//read current particle angles

		LAMMPS_NS::AtomVecEllipsoid *ellip = (LAMMPS_NS::AtomVecEllipsoid*)wrapper.lmp->atom->avec;
		TriangleSet<T> OutputMesh;

		/* temperature lattice collision and stream, apply BC iterations */
		Temperature_lattice->executeInternalProcessors();// Execute all processors and communicate appropriately


		for (plint j=0; j < nparticles; j++) {
			/* read particle positions */
			T theta = acos(ellip -> bonus[wrapper.lmp->atom->ellipsoid[j]].quat[0])*2;
			Array<T,3> normedAxis = {ellip -> bonus[wrapper.lmp->atom->ellipsoid[j]].quat[1]/sin(theta/2),
                                         	 ellip -> bonus[wrapper.lmp->atom->ellipsoid[j]].quat[2]/sin(theta/2),
                                        	 ellip -> bonus[wrapper.lmp->atom->ellipsoid[j]].quat[3]/sin(theta/2)};
			Array<T,3> position={wrapper.lmp->atom->x[j][0],wrapper.lmp->atom->x[j][1],wrapper.lmp->atom->x[j][2]};


			TriangleSet<T> SphereMesh = TriangleSet<T>("symmetricSphereSurface.stl", (T) 1e-7);

			SphereMesh.scale(r[0],r[1],r[2]);
			SphereMesh.rotateAtOrigin(normedAxis, theta);
			SphereMesh.translate(position);
			
			OutputMesh.append(SphereMesh);

			SphereMesh.scale(units.getLbLength(1),units.getLbLength(1),units.getLbLength(1));
			//OblateMesh.translate({0.5,0.5,0.5});//parameters.getNx() built an extra layer for the domain, which makes the translate(position) not exact at center.
			DEFscaledMesh<T> *MeshDef = new DEFscaledMesh<T>(SphereMesh, 0, 0, 0, Dot3D(0, 0, 0));


			std::vector<Array<T,3> > vertices;
			std::vector<Array<T,3> > VertexNormal;
			std::vector<T> areas; 
			std::vector<Array<T,3> > OuterVertices;
			std::vector<Array<T,3> > InnerVertices;
			std::vector<T> HeatFlux(heatflux.NumberOfElements());

			/* calculate surface heat fluxes */
			for (pluint iVertex = 0; iVertex < (pluint) MeshDef->getMesh().getNumVertices(); iVertex++) {
				vertices.push_back(MeshDef->getMesh().getVertex(iVertex));
				areas.push_back(MeshDef->getMesh().computeVertexArea(iVertex));
				VertexNormal.push_back(MeshDef->getMesh().computeVertexNormal(iVertex, true));
					
				OuterVertices.push_back(vertices.back() + 2.25*VertexNormal.back()); 
				//effective thickness 1.25*DeltaX + 1*DeltaX for the first-order one-sided difference approximations (Suzuki2016)
				InnerVertices.push_back(vertices.back() + 1.25*VertexNormal.back());
				
				if(OuterVertices.back()[0] > (parameters.getNx()-1)) {
					OuterVertices.back()[0] -= (parameters.getNx()-1);
				}
				if(OuterVertices.back()[0] < 0) {
					OuterVertices.back()[0] += (parameters.getNx()-1);
				}

				if(InnerVertices.back()[0] > (parameters.getNx()-1)) {
					InnerVertices.back()[0] -= (parameters.getNx()-1);
				}
				if(InnerVertices.back()[0] < 0) {
					InnerVertices.back()[0] += (parameters.getNx()-1);
				}

//if(j=40) {
//pcout << "OuterVertices: " << OuterVertices.back()[0] << "," << OuterVertices.back()[1] << "," << OuterVertices.back()[2] << std::endl;
//}
				T OuterVerticesTemperature;
				T InnerVerticesTemperature;
				if(OuterVertices.back()[1] > (parameters.getNy()-1) || OuterVertices.back()[1] < 0 || OuterVertices.back()[2] > (parameters.getNz()-1) || OuterVertices.back()[2] < 0) {
					OuterVerticesTemperature = T_wall;
				} else {
					OuterVerticesTemperature = TrilinearInterpolation(Temperature_lattice, OuterVertices.back());
				}

				if(InnerVertices.back()[1] > (parameters.getNy()-1) || InnerVertices.back()[1] < 0 || InnerVertices.back()[2] > (parameters.getNz()-1) || InnerVertices.back()[2] < 0) {
					InnerVerticesTemperature = T_wall;
				} else {
					InnerVerticesTemperature = TrilinearInterpolation(Temperature_lattice, InnerVertices.back());
				}

				HeatFlux[iVertex] = (OuterVerticesTemperature-InnerVerticesTemperature) / units.getPhysLength(1) * k_f; //positive value if with the direction inwards
			}
			delete MeshDef;

			heatflux.SetData(&HeatFlux[0],heatflux.NumberOfElements());
			
			VertexNormal.clear();
			OuterVertices.clear();
			InnerVertices.clear();
			HeatFlux.clear();

			/* read particle temperature */
			ParticleTemperatures.SetData(&ParticleTemperature[j*4357],ParticleTemperatures.NumberOfElements());

			/* run matlab */
			IT(1,1) = iT;

			particleThermalComputation(2, SurfaceTemperatures, ParticleTemperatures, IT, VtkSteps, Dt_phys, heatflux, ParticleTemperatures);//2 output; 5 inputs.
			
			/* save particle temperature */
			ParticleTemperatures.GetData(&ParticleTemperature[j*4357],ParticleTemperatures.NumberOfElements());

			// Instantiate the immersed wall data and performed the immersed boundary iterations.
			instantiateImmersedWallData(vertices, areas, container);
			vertices.clear();
			areas.clear(); 
			for (int i = 0; i < ibIter; i++) {
				indexedInamuroAdvectionDiffusionIteration(SurfaceTemperatureFunc(SurfaceTemperatures),
				*temperatureBar, container, (T) 1.0 / omega_temperature);
			}
pcout << "computing particles: " << j+1 << "/" << nparticles << std::endl;

		}


		
		if(iT%vtkSteps == 0 && iT > 0) {
			pcout << "Creating the immersed particle surfaces at " << iT  <<  std::endl;
			OutputMesh.writeBinarySTL(createFileName(demOutDir+"surface_", iT, 6) + ".stl");
		}



		/*if(iT%logSteps == 0){
			end = clock();
			T time = difftime(end,loop)/((T)CLOCKS_PER_SEC);
			T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
			T mlups = ((T) (parameters.getNx()*parameters.getNy()*parameters.getNz()*logSteps))/time/1e6;
			pcout << "time: " << time << " " ;
			pcout << "calculating at " << mlups << " MLU/s"
			      << " | total time running: " << totaltime << std::endl;
			loop = clock();
		} */

	}


	// terminate the lib
	libMatlabFETerminate();     
	// terminate MCR  
	mclTerminateApplication();

	delete j;
	delete temperatureBar;
	delete lattice;
	delete Temperature_lattice;


	T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
	T totalmlups = ((T) (parameters.getNx()*parameters.getNy()*parameters.getNz()*(maxSteps+1)))/totaltime/1e6;
	pcout << " ********************** \n"
	      << "total time: " << totaltime
	      << " calculating at " << totalmlups << " MLU/s" << std::endl;
}
