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
				TBTKExit(
					"TCallBack::getHoppingAmplitude()",
					"Should never happened",
					""
				);
			} 
		}
		else{
			const Geometry& geometry = model.getGeometry();
			Vector3d toCoordinate = geometry.getCoordinate(to);
			Vector3d fromCoordinate = geometry.getCoordinate(from);
			Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);

			Vector3d k({
				((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2))),
				((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2))),
				0
			});
			return -t*f(to, from)*exp(i*Vector3d::dotProduct(k, distanceMinimizingTranslation));
			TBTKExit(
				"TCallBack::getHoppingAmplitude()",
				"Should never happened",
				""
			);
		} 
	}

	Vector3d getDistanceMinimizingTranslation(const Vector3d &toCoordinate, const Vector3d &fromCoordinate) const {
		double smallestNorm = std::numeric_limits<double>::infinity(); 		
		Vector3d requiredTranslation({0,0,0});
		for (int x=-1; x<2; x++){
			for (int y=-1; y<2; y++){			
				Vector3d translation = x*unitCellBasis[0]+ y*unitCellBasis[1];
				Vector3d translatedFromCoordinate = fromCoordinate + translation;
				Vector3d difference = toCoordinate - translatedFromCoordinate;

				double norm = difference.norm();
	
				if(smallestNorm > norm){
					smallestNorm = norm;

					requiredTranslation = translation;
				}
			};
		};

		return requiredTranslation;
	}

	// This function is responsable to handle with one pair of indices
	complex<double> f(const Index &to, const Index &from) const{
		Vector3d toCoordinate=model.getGeometry().getCoordinate(to);
		Vector3d fromCoordinate=model.getGeometry().getCoordinate(from);

		Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);
		Vector3d difference = fromCoordinate + distanceMinimizingTranslation - toCoordinate;
		double smallestDistance = difference.norm();

		double beta=3; 		

		if (to[2] == from[2]){
			return exp(-beta*(smallestDistance/(a/sqrt(3)) - 1));
		}
		else{
			return Vector3d::dotProduct(difference.unit(), {0, 0, 1})*exp(-beta*(smallestDistance/(a/sqrt(3)) - 1));
		}
	}

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

class Graphene{
public:
	Graphene(unsigned int SIZE_K):
	//Reciprocal-Space - BrillouinZone
	SIZE_K(SIZE_K),
	SIZE_KX(SIZE_K),
	SIZE_KY(SIZE_K),
	K_POINTS_PER_PATH(SIZE_K/2),
 	//RealSpace - LatticeInformation
	a(2.45),
	numAtomsUnitCell(4),
	SIZE_X({5,5}),
	SIZE_Y({5,5}),
	abStacking(true),
	rectangularUnitCell(true),
	unitCellSize({SIZE_X[0]*a*sqrt(3), SIZE_Y[0]*a}),
	unitCellBasis({{unitCellSize[0], 0, 0},{0, unitCellSize[1], 0}}),
	//Hamiltonian Parameters.
	t(1.0),
	tCallBack(a, model),
	//Setup energy window
	LOWER_BOUND(-2),
	UPPER_BOUND(+2),
	RESOLUTION(1000),
	//Gaussian smoothing parameters
	SMOOTHING_SIGMA(0.1),
	SMOOTHING_WINDOW(51)
	{
	}

//Reciprocal-Space - BrillouinZone
	vector<vector<int>> generateAllPoints(){
		vector<vector<int>> allKPoints;
		for (int kx=0; kx<SIZE_KX; kx++){
			for (int ky=0; ky<SIZE_KY; ky++){
				allKPoints.push_back(vector<int>());
				allKPoints.back().push_back(kx);
				allKPoints.back().push_back(ky);
			}
		}
		return allKPoints;
	}
	vector<vector<vector<int>>> generateKPaths(){
		vector<vector<vector<int>>> paths;
		for(unsigned int p = 0; p < 5; p++){
			paths.push_back(vector<vector<int>>());
			vector<vector<int>> &path = paths.back();
			//Loop over a single path.
			for(unsigned int n = 0; n < (unsigned int)K_POINTS_PER_PATH; n++){
				
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
					kx = K_POINTS_PER_PATH;
					ky = n;
					break;
				
				//M-Y
				case 2: 
					kx = K_POINTS_PER_PATH-n;
					ky = K_POINTS_PER_PATH;
					break;
				//Y-G
				case 3: 
					kx = 0;
					ky = (K_POINTS_PER_PATH-n);
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
	vector<vector<int>> kPathsToKpoints(const vector<vector<vector<int>>> &kPaths){
		vector<vector<int>> kPoints;
		for(auto path : kPaths){
			kPoints.insert(kPoints.end(), path.begin(), path.end()); 
		}

		// remove repeated numbers
		for(unsigned int n=0; n<kPoints.size(); n++){
			int kx = kPoints[n][0];
			int ky = kPoints[n][1];

			for(unsigned int c=n+1; c<kPoints.size(); c++){

				if(kx==kPoints[c][0] && ky==kPoints[c][1]){
					kPoints.erase(kPoints.begin()+c, kPoints.begin()+c+1);
					c--;
				}
			}
		}

		return kPoints;
	}

//Model
	void setupModel(vector<vector<int>> &kPoints){	
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

		addHoppingAmplitude(kPoints);

		model.construct();
		Timer::tock();
	}	
	void addHoppingAmplitude(const vector<vector<int>> &kPoints){
		
		for (vector<int> k:kPoints){
			int kx=k[0];
			int ky=k[1];

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

			for(int toX = 0; toX < SIZE_X[0]; toX++){
				for (int toY = 0; toY < SIZE_Y[0]; toY++){					
					for (int toSite=0; toSite<4; toSite++){
						for(int fromX = 0; fromX < SIZE_X[1]; fromX++){
							for (int fromY = 0; fromY < SIZE_Y[1]; fromY++){					
								for (int fromSite=0; fromSite<4; fromSite++){
										model << HoppingAmplitude(tCallBack, {kx, ky, 0, toX, toY, toSite}, {kx, ky, 1, fromX, fromY,fromSite});
								}
							}
						}
					}
				}
			} 
		}
	}

//RealSpace - LatticeInformation
	void setupGeometry(){
		if (rectangularUnitCell==true){
			setupRectangularGeometry();
		} else{
			setupHexagonalGeometry();
		};
	}
	void setupRectangularGeometry(){
		Timer::tick("Setup geometry");
		vector<vector<double>> deformation = {{1, 1},{SIZE_X[0]/((double)SIZE_X[1]), SIZE_Y[0]/((double)SIZE_Y[1])}};
		double ax[2]={3*(a/sqrt(3)), 3*(a/sqrt(3))};
		double ay[2]={a,a};
		double carbonLayerSeparation = 2.7;
		Vector3d R_layer({0,0, carbonLayerSeparation});
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

		if (abStacking==true){
			Vector3d layerTranslation[2]={
				Vector3d({ 0, 0, 0}),
				Vector3d({ap, 0, 0})
			};

			for (int layer=0; layer<2; layer++){
				for (int site=0; site<numAtomsUnitCell; site++){
					offSetSite[layer][site] = offSetSite[layer][site] + layerTranslation[layer];
				}
			}
		}

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
	void setupHexagonalGeometry(){
		
		TBTKExit(
			"SetupHexagonalGeometry()",
			"Still building the function",
			"");
		exit(1);
	}
	void printGeometry(){
		Timer::tick("Print geometry");
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

//Run the Model
	void setupAndRunSolver(){
		//Setup and run Solver
		Timer::tick("Run the solver");
		// Solver::BlockDiagonalizer solver;
		solver.setModel(model);
		solver.run();
		Timer::tock();
	}

//BandStructure
	void runBandStructureCalculation(){
		vector<vector<int>> kPoints = kPathsToKpoints(generateKPaths());
		numKpoints = kPoints.size();
		setupModel(kPoints);
		setupGeometry();
		printGeometry();
		setupAndRunSolver();
		calculateBandStructure();
	}
	void calculateBandStructure(){
		//BandStructure
		Timer::tick("Calculate BandStructure");
		//Property Extractor
	
		PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	
		unsigned int numBands=1*(model.getBasisSize()/numKpoints);
		Array<double> bandStructure({numBands, ((unsigned int) 5)*K_POINTS_PER_PATH}, 0);
		vector<vector<vector<int>>> paths=generateKPaths();

		unsigned int bandPlotCounter=0;
		for(auto path:paths){
			for(auto k:path){
				for (unsigned int band=0; band<numBands; band++){
					bandStructure[{band, bandPlotCounter}] = propertyExtractor.getEigenValue({k[0], k[1]}, band);
				}
				bandPlotCounter++;
			}
		}
	
		double min = bandStructure[{0,0}];
		double max = bandStructure[{numBands-1,0}];

		for(unsigned int n=0; n<((unsigned int) 5)*K_POINTS_PER_PATH; n++){
			if(min > bandStructure[{0, n}])
				min = bandStructure[{0, n}];
			if(max < bandStructure[{numBands-1,n}])
				max = bandStructure[{numBands-1,n}];
		}	
		
		//Plotting
		Plotter plotter;		
		plotter.clear();
		plotter.setLabelX("k");
		plotter.setLabelY("Energy");
		plotter.setBoundsY(-1,1);
		for (int band=0; band<(int)numBands; band++){
			plotter.plot(bandStructure.getSlice({band, _a_}), {{"color", "black"}, {"linestyle", "-"}});
		}
		for (unsigned int n=0; n<5; n++){
			plotter.plot({((double) n)*K_POINTS_PER_PATH, ((double) n)*K_POINTS_PER_PATH},{min, max},{{"color", "black"},{"linestyle","-"}});
		}
		plotter.save("figures/bandStructure.png");
		plotter.clear();
		Timer::tock();
	}

//DensityOfState
	void runDOSCalculation(){
			vector<vector<int>> kPoints = generateAllPoints();
			setupModel(kPoints);
			setupGeometry();
			printGeometry();
			setupAndRunSolver();
			calculateDOS();
		}
	void calculateDOS(){
		Timer::tick("Calculate DOS");
		
		//Property Extractor
		PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
		//Setting the energyWindow
		propertyExtractor.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);
		Property::DOS dos = propertyExtractor.calculateDOS();
		dos = Smooth::gaussian(dos, SMOOTHING_SIGMA, SMOOTHING_WINDOW);
		//Plotting
		Plotter plotter;
		plotter.plot(dos);
		plotter.save("figures/DOS.png");
		plotter.clear();
		Timer::tock();
	}

private:
	//Reciprocal-Space - BrillouinZone
	const int SIZE_K;
	const int SIZE_KX;
	const int SIZE_KY;
	const int K_POINTS_PER_PATH;
	unsigned int numKpoints;
 	//Initializing the model
	Model model;
	Solver::BlockDiagonalizer solver;
 	//RealSpace - LatticeInformation
	double a; 
 	const int numAtomsUnitCell;
	vector<int> SIZE_X;
	vector<int> SIZE_Y;
	bool abStacking;
	bool rectangularUnitCell;
	vector<double> unitCellSize;
	vector<Vector3d> unitCellBasis;
	//Hamiltonian Parameters.
	complex<double> t;
	TCallBack tCallBack;
	//Setup energy window
	const double LOWER_BOUND;
	const double UPPER_BOUND;
	const int RESOLUTION;
	//Gaussian smoothing parameters
	const double SMOOTHING_SIGMA;
	const unsigned int SMOOTHING_WINDOW;
};

int main(int argc, char **argv){
	///Initialize TBTK
	Initialize();

	//Units
	UnitHandler::setScales({"1 rad", "1 C","1 pcs","1 eV","1 Ao","1 K", "1 s"});

	//Initialize Graphene
	// Graphene grapheneDOS(12);
	// grapheneDOS.runDOSCalculation();

	Graphene grapheneBandStructure(36);
	grapheneBandStructure.runBandStructureCalculation();

	return 0;
};