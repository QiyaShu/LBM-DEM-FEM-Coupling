/*
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

/*
 * This example consists of a cubic domain that is filled with 
 * a single spherical particle. Periodic boundary conditions are applied. 
 * The particles rotates in the center under shear.
 * Heat flux by a certain temperature difference is computed.
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
#define DESCRIPTOR descriptors::D3Q19Descriptor
#define BASEDYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
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
	vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", units.getPhysVel(1));
	vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1.0/units.getPhysTime(1));

	std::unique_ptr<MultiScalarField3D<T>> temperature = computeDensity(Temperature_lattice);
	vtkOut.writeData<float>(*temperature, "temperature", 1.0);

	/*MultiScalarField3D<T> tmp(lattice);
	IBscalarQuantity sf = SolidFraction;
	applyProcessingFunctional(new GetScalarQuantityFromDynamicsFunctional<T,DESCRIPTOR,T>(sf),
                            lattice.getBoundingBox(),lattice,tmp);
	vtkOut.writeData<float>(tmp,"solidfraction",1. ); */
}


template<class TemperatureBlockLatticeT>
void writetopbotheatflux (TemperatureBlockLatticeT& Temperature_lattice,
		      T k_f, plint iter, PhysUnits3D<T> const& units,
		      const char *topfilename, const char *bottomfilename){
	std::vector<T> topheatGradient;
	for (plint dx=0; dx <= Temperature_lattice->getBoundingBox().x1-1; dx++){
		for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-1; dz ++){
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
		for(plint dz =0; dz <= Temperature_lattice->getBoundingBox().z1-1; dz ++){
			Cell<T,TEMPERATURE_DESCRIPTOR> const& cell_sample1 = Temperature_lattice->get(dx,Temperature_lattice->getBoundingBox().y0 + 1,dz);
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

	pcout<<"writing the heat flux data" << std::endl;
 }


template<class TemperatureBlockLatticeT>
T TrilinearInterpolation(TemperatureBlockLatticeT& TempLattice, Array<T,3> OuterVertices){
	T OuterVerticesTemperature;
	
	T dx = OuterVertices[0]-(int)OuterVertices[0];
	T dy = OuterVertices[1]-(int)OuterVertices[1];
	T dz = OuterVertices[2]-(int)OuterVertices[2];

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
	temperature = t.Get(globalVertexIds,1);

        return temperature;
    }
private:
    mwArray t;
};



class iniVelocityFunc {
public:
    iniVelocityFunc (T rho_, T u_, IncomprFlowParam<T> parameters_)
        : rho_f(rho_), umax(u_), parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T& rho, Array<T,3>& u) const {
        rho = rho_f;
	u[0] = umax*(((T)iY-parameters.getNy()/2)/parameters.getNy()*2);
	u[1] = 0;
	u[2] = 0;
    }
private:
    T rho_f, umax;
    IncomprFlowParam<T> parameters;
};



class iniTemperatureFunc {
public:
    iniTemperatureFunc (T tempTop_, T tempBottom_, IncomprFlowParam<T> parameters_)
        : tempTop(tempTop_), tempBottom(tempBottom_), parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T& rho, Array<T,3>& u) const {
        rho = (tempTop-tempBottom)*((T)iY/parameters.getNy());
	u[0] = 0;
	u[1] = 0;
	u[2] = 0;
    }
private:
    T tempTop, tempBottom;
    IncomprFlowParam<T> parameters;
};



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

	T L = 2e-4; //domain length
	const T lx = L, ly = L, lz = L;
	T d0 = (5.11328272e-5)*2; //particle diameter
	T N = 30; //resolution

	T shear = 1500; //shear rate
	T u0 = d0*shear;

	T dx = 1.0/N;
	T dt = 1.0/1000; // time resolution
	T u0_LB = dt/dx; 
	
	T t_0 = d0/u0;//1.0/shear; // reference time
	T dt_phys = t_0 * dt;
        T t_total = 1; // total simulated physical time space

	// this is equivalent to the variable command in LIGGGHTS/LAMMPS
	double r[3] = {d0/2, d0/2, d0/2};
	double k[3] = {k_s, k_s, k_s};//thermal conduction in different directions

	PhysUnits3D<T> units(d0,u0,nu_f,lx,ly,lz,N,u0_LB,rho_f);

	IncomprFlowParam<T> parameters(units.getLbParam());


	const plint maxSteps = ceil(t_total/dt_phys);
	const plint vtkSteps = 1000;
 	int ibIter = 5; // Iterations for the immersed boundary method.	

	//temperatur parameters setting 
	T Prandtl = nu_f/(k_f/(rho_f*cp_f));
	T LatticeDiff_O = parameters.getLatticeNu() /Prandtl;
	T LatticeDiff = (nu_f*units.getLbLength(1)*units.getLbLength(1)/units.getLbTime(1)) /Prandtl;

	T omega_temperature = (T) 1 / (TEMPERATURE_DESCRIPTOR<T>::invCs2 * LatticeDiff + (T) 0.5); //Relaxation parameter for the temperature
	
	pcout << " ---------------------------------------------- \n"
	      << "omega: " << parameters.getOmega() << "\n"
	      << "omega_temperature: " << omega_temperature << "\n"
	      << "dt_phys: " << dt_phys << "\n"
	      << "shear rate: " << shear << "\n"
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
	wrapper.setVariable("center",L/2);
	wrapper.setVariable("dx_part",r[0]*2);
	wrapper.setVariable("dy_part",r[1]*2);
	wrapper.setVariable("dz_part",r[2]*2);
	wrapper.setVariable("mass",rho_s*(4*M_PI*r[0]*r[1]*r[2]/3));
	//wrapper.setVariable("density",rho_s);

	wrapper.execFile("in.lbdem");


/*
 * Set up velocity lattice and temperature lattice for fluid
 */
	T omega = parameters.getOmega();

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
	

    	// temperature lattice 
	MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR > *Temperature_lattice = new MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR>(
		parameters.getNx(),parameters.getNy(),parameters.getNz(), 
		new AdvectionDiffusionRLBdynamics<T, TEMPERATURE_DESCRIPTOR>(omega_temperature) );
	Temperature_lattice->toggleInternalStatistics(false);
    
	MultiScalarField3D<T> *temperatureBar = generateMultiScalarField<T>((MultiBlock3D&) *Temperature_lattice, largeEnvelopeWidth).release();
	temperatureBar->toggleInternalStatistics(false);
   
	MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *Temperature_lattice, largeEnvelopeWidth).release();
	j->toggleInternalStatistics(false);

        

	// Velocity boundary conditions.
	pcout << "Generating boundary conditions." << std::endl;
	OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
            = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

	const plint nx = parameters.getNx();
	const plint ny = parameters.getNy();
	const plint nz = parameters.getNz();
	Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
	Box3D bottomLid = Box3D(0, nx-1, 0, 0, 0, nz-1);
	Box3D everythingButLid = Box3D(0, nx-1, 1, ny-2, 0, nz-1);

	boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, bottomLid,   boundary::dirichlet);
	boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, topLid,      boundary::dirichlet);

	T u = units.getLbVel(shear*L/2);//parameters.getLatticeU()*parameters.getNy()/parameters.getResolution()/2; 
	//initializeAtEquilibrium(*lattice, everythingButLid, units.getLbDensity(rho_f), Array<T,3>((T)0.,(T)0.,(T)0.) );
	//initializeAtEquilibrium(*lattice, topLid, units.getLbDensity(rho_f), Array<T,3>(u,(T)0.,(T)0.) );
	//initializeAtEquilibrium(*lattice, bottomLid, units.getLbDensity(rho_f), Array<T,3>(-u,(T)0.,(T)0.) );
	initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), iniVelocityFunc(units.getLbDensity(rho_f),u,parameters));
	setBoundaryVelocity(*lattice, topLid, Array<T,3>(u,(T)0.,(T)0.) );
	setBoundaryVelocity(*lattice, bottomLid, Array<T,3>(-u,(T)0.,(T)0.) );
	lattice->periodicity().toggle(0,true);
	lattice->periodicity().toggle(2,true);// (0 is for X, 1 for Y, 2 for Z).

	lattice->initialize(); 
	delete boundaryCondition;


	// Temperature boundary conditions
	plint temperatureTimeFactor = 1;// Ratio between the temperature and the fluid time steps
	T T_top = 1;//topLid temperature
	T T_bottom = 0;//bottomLid temperature
	T T_everythingButLid = 0.5;//flow temperature

	OnLatticeAdvectionDiffusionBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR> *bc = createLocalAdvectionDiffusionBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR>();
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, topLid, boundary::dirichlet);
	bc-> setTemperatureConditionOnBlockBoundaries(*Temperature_lattice, bottomLid, boundary::dirichlet);
    
	setBoundaryDensity(*Temperature_lattice, bottomLid, T_bottom);
	setBoundaryDensity(*Temperature_lattice, topLid, T_top);
	Temperature_lattice->periodicity().toggle(0,true);
	Temperature_lattice->periodicity().toggle(1,false);
	Temperature_lattice->periodicity().toggle(2,true);

        temperatureBar->periodicity().toggle(0, true);
        temperatureBar->periodicity().toggle(1, false);
        temperatureBar->periodicity().toggle(2, true);

        j->periodicity().toggle(0, true);
        j->periodicity().toggle(1, false);
        j->periodicity().toggle(2, true);

	//initializeAtEquilibrium(*Temperature_lattice, everythingButLid, T_everythingButLid, Array<T,3>((T)0,(T)0,(T)0));
	initializeAtEquilibrium(*Temperature_lattice, topLid, T_top, Array<T,3>((T)0,(T)0,(T)0));
	initializeAtEquilibrium(*Temperature_lattice, bottomLid, T_bottom, Array<T,3>((T)0,(T)0,(T)0));
	initializeAtEquilibrium(*Temperature_lattice, everythingButLid, iniTemperatureFunc(T_top,T_bottom,parameters));
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
		std::cout << "Could not initialize libMatlabFE!" << std::endl;  
		return -1;  
	}

	//allocate space for variables for matlab
	mwArray SurfaceTemperatures(658, 1, mxDOUBLE_CLASS);
	mwArray ParticleTemperatures(1253, 1, mxDOUBLE_CLASS);
	mwArray heatflux(658, 1, mxDOUBLE_CLASS);// perpendicular heat flux of each node on the surface
	//mwArray Radius(3, 1, mxDOUBLE_CLASS);
	//mwArray Rho(1, 1, mxDOUBLE_CLASS);
	//mwArray Cp(1, 1, mxDOUBLE_CLASS);
	//mwArray K(3, 1, mxDOUBLE_CLASS);
	mwArray IT(1, 1, mxINT32_CLASS);
	mwArray VtkSteps(1, 1, mxINT32_CLASS);
	mwArray Dt_phys(1, 1, mxDOUBLE_CLASS);

	//Radius.SetData(r,3);
	//Rho(1,1) = rho_s;
	//Cp(1,1) = cp_s;
	//K.SetData(k,3);
	VtkSteps(1,1) = vtkSteps;
	Dt_phys(1,1) = dt_phys;

	std::vector<T> initialT(ParticleTemperatures.NumberOfElements(),0.5);
	ParticleTemperatures.SetData(&initialT[0],ParticleTemperatures.NumberOfElements()); //set initial temperature for the paricle
	std::vector<T> initialSurfaceT(SurfaceTemperatures.NumberOfElements(),0.5);
	SurfaceTemperatures.SetData(&initialSurfaceT[0],SurfaceTemperatures.NumberOfElements()); //set initial temperature for the surface of particle
	//particleThermalConductionSetup(1, ParticleTemperatures, Radius, Rho, Cp, K);//1 output: ParticleTemperature; 4 input: Radius, Rho, Cp, K.


	TriangleSet<T> SphereMesh = TriangleSet<T>("sphereSurface.stl", (T) 1e-6);

	ofstream topfluxinfile;
	topfluxinfile.open("averagetopflux.txt",ios::trunc);
	topfluxinfile.close();
	ofstream bottomfluxoutfile;
	bottomfluxoutfile.open("averagebottomflux.txt",ios::trunc);
	bottomfluxoutfile.close();


/*
 * Main loops
 */
	pcout << "Staring simulation."<<std::endl;
	for (plint iT=0; iT<=maxSteps; ++iT) {

		//static bool initWithVel = true;
		//setSpheresOnLattice(*lattice,wrapper,units,false);
		setEllipsoidOnLattice(*lattice,wrapper,units,false);
		//if(initWithVel) initWithVel = false;

		if(iT%vtkSteps == 0 && iT > 0) { // LIGGGHTS does not write at timestep 0
			writeVTK(*lattice, *Temperature_lattice, parameters, units, iT);
			//SphereMesh.writeBinarySTL(createFileName(demOutDir+"surface_", iT, 6) + ".stl");
			writetopbotheatflux (Temperature_lattice, k_f, iT, units, "averagetopflux.txt", "averagebottomflux.txt");
		}

		/* Velocity lattice collision and stream, run liggghts */
		lattice->collideAndStream();


		getForcesFromLattice(*lattice,wrapper,units);
  
		// equilvalent to the "run" command in LIGGGHTS/LAMMPS
		wrapper.run(demSubsteps);

		//read current particle angles

		LAMMPS_NS::AtomVecEllipsoid *ellip = (LAMMPS_NS::AtomVecEllipsoid*)wrapper.lmp->atom->avec;

		T theta = acos(ellip -> bonus[0].quat[0])*2;
		Array<T,3> normedAxis = {ellip -> bonus[0].quat[1]/sin(theta/2),
                                         ellip -> bonus[0].quat[2]/sin(theta/2),
                                         ellip -> bonus[0].quat[3]/sin(theta/2)};

		Array<T,3> position={wrapper.lmp->atom->x[0][0],wrapper.lmp->atom->x[0][1],wrapper.lmp->atom->x[0][2]};

		TriangleSet<T> SphereMesh = TriangleSet<T>("sphereSurface.stl", (T) 1e-6);
		std::vector<Array<T,3> > vertices;
		std::vector<Array<T,3> > VertexNormal;
		std::vector<T> areas; 
		std::vector<Array<T,3> > OuterVertices;
		std::vector<T> HeatFlux(heatflux.NumberOfElements());

		SphereMesh.scale(r[0],r[1],r[2]);
		SphereMesh.rotateAtOrigin(normedAxis, theta);
		SphereMesh.translate(position);

		if(iT%vtkSteps == 0 && iT > 0) {
			SphereMesh.writeBinarySTL(createFileName(demOutDir+"surface_", iT, 6) + ".stl");
			pcout << "Creating the immersed particle surface at " << iT << " steps." << std::endl;
		}


		/* temperature lattice collision and stream, apply BC iterations */

		Temperature_lattice->executeInternalProcessors();// Execute all processors and communicate appropriately

		SphereMesh.scale(units.getLbLength(1),units.getLbLength(1),units.getLbLength(1));
		DEFscaledMesh<T> *MeshDef = new DEFscaledMesh<T>(SphereMesh, 0, 0, 0, Dot3D(0, 0, 0));
		//pcout << "The particle surface has " << MeshDef->getMesh().getNumVertices() << " vertices and " <<
		//	MeshDef->getMesh().getNumTriangles() << " triangles." << std::endl;


		std::vector<T> SurfaceTemperature(SurfaceTemperatures.NumberOfElements());
		SurfaceTemperatures.GetData(&SurfaceTemperature[0],SurfaceTemperatures.NumberOfElements());

		for (pluint iVertex = 0; iVertex < (pluint) MeshDef->getMesh().getNumVertices(); iVertex++) {
			vertices.push_back(MeshDef->getMesh().getVertex(iVertex));
			areas.push_back(MeshDef->getMesh().computeVertexArea(iVertex));
			VertexNormal.push_back(MeshDef->getMesh().computeVertexNormal(iVertex, true));
						
			OuterVertices.push_back(vertices.back() + 2.25*VertexNormal.back()); 
			//effective thickness 1.25*DeltaX + 1*DeltaX for the first-order one-sided difference approximations (Suzuki2016)

			T OuterVerticesTemperature = TrilinearInterpolation(Temperature_lattice, OuterVertices.back());
			HeatFlux[iVertex] = (OuterVerticesTemperature-SurfaceTemperature[iVertex]) / units.getPhysLength(1) * k_f; //positive value if with the direction inwards		
		}
		delete MeshDef;


		/* after reading heat flux, run MATLAB  */

		IT(1,1) = iT;
		heatflux.SetData(&HeatFlux[0],heatflux.NumberOfElements());

		particleThermalComputation(2, SurfaceTemperatures, ParticleTemperatures, IT, VtkSteps, Dt_phys, heatflux, ParticleTemperatures);//2 output; 5 inputs.


		// Instantiate the immersed wall data and performed the immersed boundary iterations.
		instantiateImmersedWallData(vertices, areas, container);
		for (int i = 0; i < ibIter; i++) {
			indexedInamuroAdvectionDiffusionIteration(SurfaceTemperatureFunc(SurfaceTemperatures),
			*temperatureBar, container, (T) 1.0 / omega_temperature);
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
	//delete lattice;


	T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
	T totalmlups = ((T) (parameters.getNx()*parameters.getNy()*parameters.getNz()*(maxSteps+1)))/totaltime/1e6;
	pcout << " ********************** \n"
	      << "total time: " << totaltime
	      << " calculating at " << totalmlups << " MLU/s" << std::endl;
}
