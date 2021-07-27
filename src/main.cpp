#include "TBTK/BrillouinZone.h"
#include "TBTK/Model.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Smooth.h"
#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const complex<double> i(0,1);

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();
	UnitHandler::setScales({"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Lattice parameters
	double a = 2.5;
	const int SIZE_X = 2;
	const int SIZE_Y = 2;

	// //Hamiltonian parameters;
	complex<double> t=1.0;
	// double gamma = 0.3; 
	// double voltage = 0.0;

	// //Auxiliar
	int SITE_0 = 0;
	int SITE_1 = 1;
	int SITE_2 = 2;
	int SITE_3 = 3;

	// int LAYER0 = 0;
	// int LAYER1 = 1;

	// int subLatticeA=0;
	// int subLatticeB=1;

	//Lattice Vector
	Vector3d latticeVector[3];
	latticeVector[0] = Vector3d({   a, 			 0,	0});
	latticeVector[1] = Vector3d({-a/2, a*sqrt(3)/2,	0});
	latticeVector[2] = Vector3d({	0,			 0, a});

	//Define NN for siteA
	Vector3d nnVector[3];
	nnVector[0] =  latticeVector[0] + 2*latticeVector[1]/3;
	nnVector[1] = -latticeVector[1] + nnVector[0];
	nnVector[2] = -latticeVector[0] - nnVector[1] + nnVector[0];

	//Calculate the reciprocal lattice vectors.
	Vector3d k[3];
	for(unsigned int n=0; n<3; n++){
		k[n] = 2*M_PI*latticeVector[(n+1)%3]*latticeVector[(n+2)%3]/(
			Vector3d::dotProduct(latticeVector[n], latticeVector[(n+1)%3]*latticeVector[(n+2)%3])
		);
	}

	//Setup BZ
	vector<unsigned int> numMeshPoints = {2*SIZE_X, 2*SIZE_Y};
	
	BrillouinZone brillouinZone({
		{k[0].x, k[0].y},{k[1].x, k[1].y}
		},
		SpacePartition::MeshType::Nodal
	);
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(numMeshPoints);
 
	//Model and hopping
	const int SIZE_KX = 2;
	const int SIZE_KY = 2;
	const double R_X  = 2*a + 2*a*0.5;
	const double R_Y  = 2*a*sqrt(3)/2;

	Model model;
	model.setVerbose(true);

	for (int kx=0; kx<SIZE_KX; kx++){
		double KX = 2*kx*M_PI/SIZE_KX;
		for (int ky=0; ky<SIZE_KY; ky++){
			double KY = 2*ky*M_PI/SIZE_KY;
			for (int x =0; x<SIZE_X; x++){
				for (int y=0; y<SIZE_Y; y++){
					for (int spin=0; spin<2; spin++){
						model << HoppingAmplitude(-t, {kx, ky, x, y, SITE_1, spin}, {kx, ky, x, y, SITE_0, spin} ) + HC;
						model << HoppingAmplitude(-t, {kx, ky, x, y, SITE_2, spin}, {kx, ky, x, y, SITE_1, spin} ) + HC;
						model << HoppingAmplitude(-t, {kx, ky, x, y, SITE_3, spin}, {kx, ky, x, y, SITE_2, spin} ) + HC;

						if(x+1 < SIZE_X){
							model << HoppingAmplitude(-t, {kx, ky, (x+1)%SIZE_X, y, SITE_1, spin}, {kx, ky, x, y, SITE_0, spin} ) + HC;
						} else {
							model << HoppingAmplitude(-t*exp(i*KX*(SIZE_X*R_X)), {kx, ky, (x+1)%SIZE_X, y, SITE_1, spin}, {kx, ky, x, y, SITE_0, spin} ) + HC;
						}

						if(y+1 < SIZE_Y){
							model << HoppingAmplitude(-t, {kx, ky, x, (y+1)%SIZE_Y, SITE_1, spin}, {kx, ky, x, y, SITE_0, spin} ) + HC;
							model << HoppingAmplitude(-t, {kx, ky, x, (y+1)%SIZE_Y, SITE_3, spin}, {kx, ky, x, y, SITE_2, spin} ) + HC;
						} else {
							model << HoppingAmplitude(-t*exp(i*KY*(SIZE_X*R_Y)), {kx, ky, x, (y+1)%SIZE_Y, SITE_1, spin}, {kx, ky, x, y, SITE_0, spin} ) + HC;
							model << HoppingAmplitude(-t*exp(i*KY*(SIZE_X*R_Y)), {kx, ky, x, (y+1)%SIZE_Y, SITE_3, spin}, {kx, ky, x, y, SITE_2, spin} ) + HC;
						}
					}
				}
			}
		}
	}

	model.construct();
	Streams::out << "1/n";

	//Setup and run Solver: Diagonalizer
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();
	Streams::out << "2/n";

	//Create PropertyExtractor
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	Streams::out << "3/n";

	/*
		ENERGY LEVEL AND DOS
	*/

	//Setup energy window
	const double LOWER_BOUND = -5;
	const double UPPER_BOUND = +5;
	const int RESOLUTION = 1000;

	propertyExtractor.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);
	Streams::out << "4\n";

	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();
	Streams::out << "5\n";

	Plotter plotter;
	plotter.plot(eigenValues);
	plotter.save("figures/eigenValues.png");
	Streams::out << "6\n";

	//DOS
	Property::DOS dos = propertyExtractor.calculateDOS();
	Streams::out << "7\n";

	const double SMOOTHING_SIGMA = 0.1;
	const unsigned int SMOOTHING_WINDOW = 51;
	dos = Smooth::gaussian(dos, SMOOTHING_SIGMA, SMOOTHING_WINDOW);
	Streams::out << "8\n";

	plotter.clear();
	plotter.plot(dos);
	plotter.save("figures/DOS.png");
	Streams::out << "9\n";

	/*
	Band Structure
	*/

	const int K_POINTS_PER_PATH = 100;
	//High Symmetry Points
	Vector3d G({     0,       0,   0});
	Vector3d M({M_PI/a, -M_PI/a,   0});
	Vector3d K({4*M_PI/(3*a),       0,   0});

	vector<vector<Vector3d>> paths = {{G,M},{M,K},{K,G}};

	Array<double> bandStructure({4, 3*K_POINTS_PER_PATH}, 0);
	Range interpolator(0,1, K_POINTS_PER_PATH);

	double min = bandStructure[{0,0}];
	double max = bandStructure[{3,0}];

	for(unsigned int n=0; n<3*K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0,n}])
			min = bandStructure[{0,n}];

		if(max < bandStructure[{3,n}])
			max = bandStructure[{3,n}];
	}

	plotter.clear();
	plotter.setLabelX("k");
	plotter.setLabelY("Energy");
	plotter.setBoundsY(-1,1);
	plotter.setBoundsX(150,250);
	plotter.plot(bandStructure.getSlice({0, _a_}), {{"color", "black"}});
	plotter.plot(bandStructure.getSlice({1, _a_}), {{"color", "black"}});
	plotter.plot(bandStructure.getSlice({2, _a_}), {{"color", "black"}});
	plotter.plot(bandStructure.getSlice({3, _a_}), {{"color", "black"}});
	plotter.save("figures/bandStructure.png");

	return 0;
}
