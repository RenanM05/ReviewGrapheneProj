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
#include <fstream>
#include <limits>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const complex<double> i(0,1);

//Setup energy window
const double LOWER_BOUND=-2;
const double UPPER_BOUND=+2;
const int RESOLUTION=1000;

//Gaussian smoothing parameters
const double SMOOTHING_SIGMA = 0.1;
const unsigned int SMOOTHING_WINDOW = 51;

//RealSpace - LatticeInformation
double a=2.5;
const int numAtomsUnitCell = 4;
vector<int> SIZE_X = {4,4};
vector<int> SIZE_Y = {4,4};
vector<double> unitCellSize = {SIZE_X[0]*a*sqrt(3), SIZE_Y[0]*a};
vector<Vector3d> unitCellBasis = {{unitCellSize[0], 0, 0},{0, unitCellSize[1], 0}};

//Reciprocal-Space - BrillouinZone
const int SIZE_K  = 12;
const int SIZE_KX = SIZE_K;
const int SIZE_KY = SIZE_K;

//Hamiltonian parameters
complex<double> t=1.0;

//Initializing the model
Model model; 
Solver::BlockDiagonalizer solver;

class TCallBack: public HoppingAmplitude::AmplitudeCallback{
public: 
	TCallBack(double a, const Model& model): a(a), model(model){
	}

	complex<double> getHoppingAmplitude(const Index &to, const Index &from) const{	
		
		int kx = to[0];
		int ky = to[1];
			
		int toLayer = to[2];
		int fromLayer = from[2];
		
		int toX = to[3];
		int fromX = from[3];
		
		int toY = to[4];
		int fromY = from[4];

		if(toLayer == fromLayer){
			int layer = toLayer;
						
			if (toX == fromX && toY == fromY){ 
				return -t*f(to, from);
			}

			// To connect plaquette
			if (toY == fromY){
				if (toX == SIZE_X[layer]-1 && fromX == 0){
					double KX = ((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2)));
					return -t*f(to, from)*exp(i*KX*(unitCellSize[0]));
				} else if (fromX == SIZE_X[layer]-1 && toX == 0){
					double KX = ((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2)));
					return -t*f(to, from)*exp(-i*KX*(unitCellSize[0]));
				} else{
					return -t*f(to, from);
				}
				
			} else if(toX == fromX){
				if (toY == SIZE_Y[layer]-1 && fromY == 0){
					double KY = ((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2)));
					return -t*f(to, from)*exp(i*KY*(unitCellSize[1]));
				} else if (fromY == SIZE_Y[layer]-1 && toY == 0){
					double KY = ((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2)));
					return -t*f(to, from)*exp(-i*KY*(unitCellSize[1]));
				} else {
					return -t*f(to, from);
				}				
			} else{
				Streams::out << "Should never happened" ;
				exit(1);
			} 
		}
		else{
			Streams::out << "Should never happened" ;
			exit(1);
		} 
	}

	complex<double> f(const Index& to, const Index& from) const{
		
		Vector3d toCoordinate=model.getGeometry().getCoordinate(to);
		Vector3d fromCoordinate=model.getGeometry().getCoordinate(from);
		
		double smallestNorm = std::numeric_limits<double>::infinity();
		for (int x=-1; x<2; x++){
			for (int y=-1; y<2; y++){
				Vector3d translatedFromCoordinate = fromCoordinate + x*unitCellBasis[0]+ y*unitCellBasis[1];
				Vector3d difference = toCoordinate - translatedFromCoordinate;

				double norm = difference.norm();
				if(smallestNorm > norm)
					smallestNorm = norm;
			}
		}
	
		double beta = 3;
		//When smallestNorm = 1.42, then exp->0. 
		//Therefore, the hopping between two carbons is -t
		return exp(-beta*(smallestNorm/(a/sqrt(3)) - 1));
	}

	// set the internal t to a newT
	void setT(complex<double> t){this->t = t;}
	void setSIZE_KX(int SIZE_KX){this-> SIZE_KX = SIZE_KX;}
	void setSIZE_KY(int SIZE_KY){this-> SIZE_KY = SIZE_KY;}
	void setSIZE_X(vector<int> SIZE_X){this-> SIZE_X = SIZE_X;}
	void setSIZE_Y(vector<int> SIZE_Y){this-> SIZE_Y = SIZE_Y;}
	void setUnitCellSize(vector<double> unitCellSize){this-> unitCellSize = unitCellSize;}
	void setUnitCellBasis(vector<Vector3d> unitCellBasis){this-> unitCellBasis = unitCellBasis;}

private:
	complex<double> t;
	vector<int> SIZE_X;
	vector<int> SIZE_Y;
	vector<double> unitCellSize;
	vector<Vector3d> unitCellBasis;
	int SIZE_KX;
	int SIZE_KY;
	double R_X;
	double R_Y;
	double a;
	const Model& model;
};

vector<vector<vector<int>>> generateKPaths(unsigned int kPointsPerPath){
	vector<vector<vector<int>>> paths;
	for(unsigned int p = 0; p < 5; p++){
		paths.push_back(vector<vector<int>>());
		vector<vector<int>> &path = paths.back();
		//Loop over a single path.
		for(unsigned int n = 0; n < kPointsPerPath; n++){
			
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
				kx = kPointsPerPath;
				ky = n;
				break;
			
			//M-Y
			case 2: 
				kx = kPointsPerPath-n;
				ky = kPointsPerPath;
				break;
			//Y-G
			case 3: 
				kx = 0;
				ky = (kPointsPerPath-n);
				break;
			//G-M
			case 4: 
				kx = n;
				ky = n;
				break;
			default:
				break;
			}
			path.push_back({kx, ky});
		}
	}
	return paths;
}

void setupModel(TCallBack &tCallBack){	
	//Setting the Model
	Timer::tick("Setup the Model");
	model.setVerbose(true);

	// TCallBack tCallBack(a, model);
	tCallBack.setT(t);
	tCallBack.setSIZE_X(SIZE_X);
	tCallBack.setSIZE_Y(SIZE_Y);
	tCallBack.setSIZE_KX(SIZE_KX);
	tCallBack.setSIZE_KY(SIZE_KY);
	tCallBack.setUnitCellSize(unitCellSize);
	tCallBack.setUnitCellBasis(unitCellBasis);

	for (int kx =0; kx < SIZE_KX; kx++){
		for (int ky =0; ky < SIZE_KY; ky++){
			for(int layer=0; layer<2; layer++){
				for(int x = 0; x < SIZE_X[layer]; x++){
					for (int y = 0; y < SIZE_Y[layer]; y++){
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, x, y, 1}, {kx, ky, layer, x,y,0}) + HC;
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, x, y, 2}, {kx, ky, layer, x,y,1}) + HC;
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, x, y, 3}, {kx, ky, layer, x,y,2}) + HC;
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, (x+1)%SIZE_X[layer],y,0}, {kx, ky, layer, x,y,3}) + HC;
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, x,(y+1)%SIZE_Y[layer],0}, {kx, ky, layer, x,y,1}) + HC;
						model << HoppingAmplitude(tCallBack, {kx, ky, layer, x,(y+1)%SIZE_Y[layer],3}, {kx, ky, layer, x,y,2}) + HC;						
					}
				}
			}
		} 
	}

	model.construct();
	Timer::tock();
}

void setupGeometry(){
	Timer::tick("Setup geometry");
	vector<vector<double>> deformation = {{1, 1},{SIZE_X[0]/((double)SIZE_X[1]), SIZE_Y[0]/((double)SIZE_Y[1])}};
	double ax[2]={3*(a/sqrt(3)), 3*(a/sqrt(3))};
	double ay[2]={a,a};
	double layerSeparation = 3.4;
	Vector3d R_layer({0,0, layerSeparation});
	double ap = a/sqrt(3);
	Vector3d offSetSite[2][4]={
		{
			Vector3d({0,0,0}),
			Vector3d({ap*0.5, ap*sqrt(3)/2, 0}),
			Vector3d({ap*0.5+ap, ap*sqrt(3)/2, 0}),
			Vector3d({2*ap*0.5+ap, 0, 0})
		 },
		{
			Vector3d({0,0,0}),
			Vector3d({ap*0.5, ap*sqrt(3)/2, 0}),
			Vector3d({ap*0.5+ap, ap*sqrt(3)/2, 0}),
			Vector3d({2*ap*0.5+ap, 0, 0})
		}
	};

	for (int layer=0; layer<2; layer++){
		for (int site=0; site<numAtomsUnitCell; site++){
			offSetSite[layer][site].x *= deformation[layer][0];
			offSetSite[layer][site].y *= deformation[layer][1];
		}
	}

	Geometry& geometry = model.getGeometry();
	for (int kx=0; kx<SIZE_KX; kx++){
		for (int ky=0; ky<SIZE_KY; ky++){
			for(int layer=0; layer<2; layer++){
				Vector3d layerPosition=layer*R_layer;
				for(int x=0; x<SIZE_X[layer]; x++){
					Vector3d xVector({x*ax[layer],0,0});
					for (int y=0; y<SIZE_Y[layer]; y++){
						Vector3d yVector({0,y*ay[layer],0});					
						for(int site=0; site<numAtomsUnitCell; site++){
							Vector3d plaquettePosition=layerPosition+deformation[layer][0]*xVector+deformation[layer][1]*yVector;
							geometry.setCoordinate({kx, ky, layer, x, y, site}, (plaquettePosition+offSetSite[layer][site]).getStdVector());
						}
					}
				}
			}
		} 
	}
	Timer::tock();
}

void printGeometry(){
	Timer::tick("Save geometry output file");
	ofstream fout0("systemInfo/geometry/coordinatesLayer0");
	ofstream fout1("systemInfo/geometry/coordinatesLayer1");
	Geometry& geometry = model.getGeometry();
	for(int layer=0; layer<2; layer++){
		for(int x=0; x<SIZE_X[layer]; x++){
			for (int y=0; y<SIZE_Y[layer]; y++){
				for(int site=0; site<numAtomsUnitCell; site++){
					Vector3d position = geometry.getCoordinate({0, 0, layer, x, y, site});
					if (layer==0)
						fout0 << position.x << "\t" << position.y << "\n"; 
					else
						fout1 << position.x << "\t" << position.y << "\n"; 
				}
			}
		}
	}
	
	fout0.close();
	fout1.close();
	Timer::tock();

}

void setupAndRunSolver(){
	//Setup and run Solver
	Timer::tick("Run the solver");
	// Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();
	Timer::tock();
}

void extractProperties(){
	Timer::tick("Extracting Properties");
	//Property Extractor
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	Timer::tock();
	
	//EnergyLevels - Eigenvalues
	propertyExtractor.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

	//DensityOfStates
	Property::DOS dos = propertyExtractor.calculateDOS();
	dos = Smooth::gaussian(dos, SMOOTHING_SIGMA, SMOOTHING_WINDOW);

	//BandStructure
	const int K_POINTS_PER_PATH = SIZE_K/2;
	Array<double> bandStructure({((unsigned int) 8)*(SIZE_X[0]*SIZE_Y[0] + SIZE_X[1]*SIZE_Y[1]), 5*K_POINTS_PER_PATH}, 0);
	vector<vector<vector<int>>> paths=generateKPaths(K_POINTS_PER_PATH);

	unsigned int bandPlotCounter=0;
	for(auto path:paths){
		for(auto k:path){
			for (unsigned int band=0; band <((unsigned int) 8)*(SIZE_X[0]*SIZE_Y[0]+SIZE_X[1]*SIZE_Y[1]); band++){
				bandStructure[{band, bandPlotCounter}] = propertyExtractor.getEigenValue({k[0], k[1]}, band);
			}
			bandPlotCounter++;
		}
	}
	
	double min = bandStructure[{0,0}];
	double max = bandStructure[{(unsigned int)8*(SIZE_X[0]*SIZE_Y[0]+SIZE_X[1]*SIZE_Y[1])-1,0}];

	for(unsigned int n=0; n<5*K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0, n}])
			min = bandStructure[{0, n}];
		if(max < bandStructure[{(unsigned int)8*(SIZE_X[0]*SIZE_Y[0]+SIZE_X[1]*SIZE_Y[1])-1,n}])
			max = bandStructure[{(unsigned int)8*(SIZE_X[0]*SIZE_Y[0]+SIZE_X[1]*SIZE_Y[1])-1,n}];
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
	for (int band=0; band<((int)8)*(SIZE_X[0]*SIZE_Y[0]+SIZE_X[1]*SIZE_Y[1]); band++){
		plotter.plot(bandStructure.getSlice({band, _a_}), {{"color", "black"}, {"linestyle", "-"}});
	}
	for (unsigned int n=0; n<5; n++){
		plotter.plot({((double) n)*K_POINTS_PER_PATH, ((double) n)*K_POINTS_PER_PATH},{min, max},{{"color", "black"},{"linestyle","-"}});
	}
	plotter.save("figures/bandStructure.png");
	plotter.clear();
}

int main(int argc, char **argv){

	///Initialize TBTK
	Initialize();

	//Units
	UnitHandler::setScales({"1 rad", "1 C","1 pcs","1 eV","1 Ao","1 K", "1 s"});

	TCallBack tCallBack(a, model);
	setupModel(tCallBack);
	setupGeometry();
	printGeometry();
	setupAndRunSolver();
	extractProperties();
	
	return 0;
};

//Plot the band structure.

/*
	Vector3d r[3];
	r[0] = Vector3d({a, 0, 0});
	r[1] = Vector3d({-a/2, a*sqrt(3)/2, 0});
	r[2] = Vector3d({0, 0, a});


	Vector3d k[3];
	for(unsigned int n=0; n<3; n++){
		k[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	};

	vector<unsigned int> numMeshPoints = {2*SIZE_X, 2*SIZE_Y};
	BrillouinZone brillouinZone({{k[0].x, k[0].y},{k[1].x, k[1].y}}, SpacePartition::MeshType::Nodal);
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(numMeshPoints);

	//High Symmetry Points
	Vector3d G({     0,       0,   0});
	Vector3d M({M_PI/a, -M_PI/a,   0});
	Vector3d K({4*M_PI/(3*a),       0,   0});
    vector<vector<Vector3d>> paths = {{G,M},{M,K},{K,G}};

	Array<double> bandStructure({num_of_states, number_of_paths*K_POINTS_PER_PATH}, 0);
	Range interpolator(0, 1, K_POINTS_PER_PATH);
*/