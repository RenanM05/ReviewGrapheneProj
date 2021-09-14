#include "TBTK/BrillouinZone.h"
#include "TBTK/Model.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Range.h"
#include "TBTK/Smooth.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/TBTK.h"
#include "TBTK/Timer.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>
#include <iostream>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const complex<double> i(0,1);

int main(int argc, char **argv){

	///Initialize TBTK
	Initialize();

	//Units
	UnitHandler::setScales({"1 rad", "1 C","1 pcs","1 eV","1 Ao","1 K", "1 s"});

	//Setup energy window
	const double LOWER_BOUND=-2;
	const double UPPER_BOUND=+2;
	const int RESOLUTION=1000;

	//Gaussian smoothing parameters
	const double SMOOTHING_SIGMA = 0.1;
	const unsigned int SMOOTHING_WINDOW = 51;

	//RealSpace - LatticeInformation
	double a=2.5;
	const int SIZE_X[2] = {5,5};
	const int SIZE_Y[2] = {5,5};
	const double R_X = 2*a + 2*a*0.5;
	const double R_Y = 2*a*sqrt(3)/2;

	//Reciprocal-Space - BrillouinZone
	const int SIZE_K  = 20;
	const int SIZE_KX = SIZE_K;
	const int SIZE_KY = SIZE_K;
	
	//Hamiltonian parameters
	complex<double> t=1.0;

	//Setting the Model
	Model model; 
	model.setVerbose(true);

	for (int kx =0; kx < SIZE_KX; kx++){
 		double KX = kx*2*M_PI/SIZE_KX;
		for (int ky =0; ky < SIZE_KY; ky++){
			double KY = ky*2*M_PI/SIZE_KY;
			for(int layer=0; layer<2; layer++){
				for(int x = 0; x < SIZE_X[layer]; x++){
					for (int y = 0; y < SIZE_Y[layer]; y++){
						for (int spin = 0; spin<2; spin++){
							model << HoppingAmplitude(-t, {kx, ky, layer, x, y, 1, spin}, {kx, ky, layer, x,y,0, spin}) + HC;
							model << HoppingAmplitude(-t, {kx, ky, layer, x, y, 2, spin}, {kx, ky, layer, x,y,1, spin}) + HC;
							model << HoppingAmplitude(-t, {kx, ky, layer, x, y, 3, spin}, {kx, ky, layer, x,y,2, spin}) + HC;
							if(x+1 < SIZE_X[layer]){
								model << HoppingAmplitude(-t, {kx, ky, layer, (x+1)%SIZE_X[layer],y,0,spin}, {kx, ky, layer, x,y,3, spin}) + HC;
							} else {
								model << HoppingAmplitude(-t*exp(i*KX*(SIZE_X[layer]*R_X)), {kx, ky, layer, (x+1)%SIZE_X[layer],y,0,spin}, {kx, ky, layer, x,y,3, spin}) + HC;
							}
							if(y+1 < SIZE_Y[layer]){
								model << HoppingAmplitude(-t, {kx, ky, layer, x,(y+1)%SIZE_Y[layer],0,spin}, {kx, ky, layer, x,y,1, spin}) + HC;
								model << HoppingAmplitude(-t, {kx, ky, layer, x,(y+1)%SIZE_Y[layer],3,spin}, {kx, ky, layer, x,y,2, spin}) + HC;
							} else {
								model << HoppingAmplitude(-t*exp(i*KY*(SIZE_Y[layer]*R_Y)), {kx, ky, layer, x,(y+1)%SIZE_Y[layer],0,spin}, {kx, ky, layer, x,y,1, spin}) + HC;
								model << HoppingAmplitude(-t*exp(i*KY*(SIZE_Y[layer]*R_Y)), {kx, ky, layer, x,(y+1)%SIZE_Y[layer],3,spin}, {kx, ky, layer, x,y,2, spin}) + HC;
							}
						}
					}
				}
			}
		} 
	}

	model.construct();

	//Setup and run Solver
	Timer::tick("Setup and run Solver");
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();
	Timer::tock();

	//Property Extractor
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

	//EnergyLevels - Eigenvalues
	propertyExtractor.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

	//DensityOfStates
	Property::DOS dos = propertyExtractor.calculateDOS();
	dos = Smooth::gaussian(dos, SMOOTHING_SIGMA, SMOOTHING_WINDOW);

	//BandStructure
	const int K_POINTS_PER_PATH = SIZE_K/4;
	Array<double> bandStructure({(unsigned int)32*SIZE_X[0]*SIZE_Y[0], 5*K_POINTS_PER_PATH}, 0);
	for(unsigned int p = 0; p < 5; p++){
		//Loop over a single path.
		for(unsigned int n = 0; n < K_POINTS_PER_PATH; n++){
			
			int kx;
			int ky;
			switch (p)
			{
			//G-X		
			case 0:
				kx = n;
				ky = 0;
				break;
			//X-M
			case 1:
				kx = SIZE_K/4;
				ky = n;
				break;
			
			//M-Y
			case 2: 
				kx = SIZE_K/4-n;
				ky = SIZE_K/4;
				break;
			//Y-G
			case 3: 
				kx = 0;
				ky = (SIZE_K/4-n);
				break;
			//G-M
			case 4: 
				kx = n;
				ky = n;
				break;
			default:
				break;
			}
			for (unsigned int band=0; band <(unsigned int)32*SIZE_X[0]*SIZE_Y[0]; band++){
				bandStructure[{band, n+p*K_POINTS_PER_PATH}]=propertyExtractor.getEigenValue({kx, ky}, band);
			}		
		}
	}

	double min = bandStructure[{0,0}];
	double max = bandStructure[{(unsigned int)32*SIZE_X[0]*SIZE_Y[0]-1,0}];

	for(unsigned int n=0; n<5*K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0, n}])
			min = bandStructure[{0, n}];
		if(max < bandStructure[{(unsigned int)32*SIZE_X[0]*SIZE_Y[0]-1,n}])
			max = bandStructure[{(unsigned int)32*SIZE_X[0]*SIZE_Y[0]-1,n}];
	}	

	//Plotting
	Plotter plotter;

	//plotter.plot(eigenValues);
	//plotter.save("figures/EigenValues.png");
	//plotter.clear();

	plotter.plot(dos);
	plotter.save("figures/DOS.png");
	plotter.clear();

	plotter.clear();
	plotter.setLabelX("k");
	plotter.setLabelY("Energy");
	plotter.setBoundsY(-1,1);
	for (int band=0; band<(unsigned int)32*SIZE_X[0]*SIZE_Y[0]; band++){
		plotter.plot(bandStructure.getSlice({band, _a_}), {{"color", "black"}, {"linestyle", "-"}});
	}
	for (unsigned int n=0; n<5; n++){
		plotter.plot({((double) n)*K_POINTS_PER_PATH, ((double) n)*K_POINTS_PER_PATH},{min, max},{{"color", "black"},{"linestyle","-"}});
	}
	plotter.save("figures/bandStructure.png");
	plotter.clear();

	return 0;
};

//Plot the band structure.

/*
	Vector3d r[3];
	r[0] = Vector3d({a, 0, 0});
	r[1] = Vector3d({-a/2, a*sqrt(3)/2, 0});
	r[2] = Vector3d({0, 0, a});

	Vector3d r_AB[3];
	r_AB[0] = (r[0]+2*r[1])/3;
	r_AB[1] = -r[1]+r_AB[0];
	r_AB[2] = -r[0]-r[1]+r_AB[0];

*/

/* 
	Vector3d k[3];
	for(unsigned int n=0; n<3; n++){
		k[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	};

	vector<unsigned int> numMeshPoints = {2*SIZE_X, 2*SIZE_Y};
	BrillouinZone brillouinZone({{k[0].x, k[0].y},{k[1].x, k[1].y}}, SpacePartition::MeshType::Nodal);
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(numMeshPoints);
*/

/*  
	//High Symmetry Points
	Vector3d G({     0,       0,   0});
	Vector3d M({M_PI/a, -M_PI/a,   0});
	Vector3d K({4*M_PI/(3*a),       0,   0});
    vector<vector<Vector3d>> paths = {{G,M},{M,K},{K,G}};

	const int K_POINTS_PER_PATH = SIZE_K/2;
	Array<double> bandStructure({4, 3*K_POINTS_PER_PATH}, 0);
	Range interpolator(0, 1, K_POINTS_PER_PATH);

	double min = bandStructure[{0,0}];
	double max = bandStructure[{3,0}];

	for(unsigned int n=0; n<3*K_POINTS_PER_PATH; n++){
	 	if(min > bandStructure[{0,n}])
	 		min = bandStructure[{0,n}];

	 	if(max < bandStructure[{3,n}])
	 		max = bandStructure[{3,n}];
	 }
*/