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

#include "Graphene.h"
#include "TCallBack.h"

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

Graphene::Graphene(unsigned int SIZE_K):
parameters(SIZE_K),
//Hamiltonian Parameters.
t(1.0),
tCallBack(parameters.a, model)
{
    grapheneUnitCell = SelectGrapheneUnitCell::Rectangular;
}

//Reciprocal-Space - BrillouinZone
//OK
vector<vector<int>> Graphene::generateAllKPoints(){
	vector<vector<int>> allKPoints; 
	for (int kx =0; kx < parameters.SIZE_KX; kx++){			
		for (int ky =0; ky < parameters.SIZE_KY; ky++){
			allKPoints.push_back(vector<int>());
			allKPoints.back().push_back(kx);
			allKPoints.back().push_back(ky);
		}
	}
	return allKPoints; 
}

//TO CHECK 
vector<vector<vector<int>>> Graphene::generateKPaths(){
    vector<vector<vector<int>>> paths;
    for(unsigned int p = 0; p < 5; p++){
        paths.push_back(vector<vector<int>>());
        vector<vector<int>> &path = paths.back();
        //Loop over a single path.
        for(unsigned int n = 0; n < (unsigned int)parameters.K_POINTS_PER_PATH; n++){

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
                kx = parameters.K_POINTS_PER_PATH;
                ky = n;
                break;

            //M-Y
            case 2: 
                kx = parameters.K_POINTS_PER_PATH-n;
                ky = parameters.K_POINTS_PER_PATH;
                break;
            //Y-G
            case 3: 
                kx = 0;
                ky = (parameters.K_POINTS_PER_PATH-n);
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
vector<vector<int>> Graphene::kPathsToKpoints(const vector<vector<vector<int>>> &kPaths){
    vector<vector<int>> kPoints;
    for(auto path : kPaths){
        kPoints.insert(kPoints.end(), path.begin(), path.end()); 
    }
    return kPoints;
}

//Model
//TO CHECK
void Graphene::setupModel(const vector<int> &kPoint){	
    //Setting the Model
    Timer::tick("Setup the Model");
    // model.setVerbose(true);

    tCallBack = TCallBack(parameters.a, model);
    // tCallBack.setT(t);
    tCallBack.setSIZE_X(parameters.SIZE_X);
    tCallBack.setSIZE_Y(parameters.SIZE_Y);
    tCallBack.setSIZE_KX(parameters.SIZE_KX);
    tCallBack.setSIZE_KY(parameters.SIZE_KY);
    tCallBack.setUnitCellSize(parameters.unitCellSize);
    tCallBack.setUnitCellBasis(parameters.unitCellBasis);

    addHoppingAmplitude(kPoint);

    model.construct();
    Timer::tock();
}	
//OK
void Graphene::addHoppingAmplitude(const vector<int> &kPoint){

    int kx=kPoint[0];
    int ky=kPoint[1];

	int startPointLayer = parameters.addGrapheneBottomLayer ? 0 : 1;
	const double EPSILON = 1e-2;
	for(int toLayer = startPointLayer; toLayer<2; toLayer++){
		for(int toX = 0; toX < parameters.SIZE_X[toLayer]; toX++){
			for (int toY = 0; toY < parameters.SIZE_Y[toLayer]; toY++){					
				for (int toSite=0; toSite<4; toSite++){				
					for(int toOrbital=3; toOrbital<4; toOrbital++){
						for(int fromLayer = startPointLayer; fromLayer<2; fromLayer++){
							for(int fromX = 0; fromX < parameters.SIZE_X[fromLayer]; fromX++){
								for (int fromY = 0; fromY < parameters.SIZE_Y[fromLayer]; fromY++){					
									for (int fromSite=0; fromSite<4; fromSite++){
										for(int fromOrbital=3; fromOrbital<4; fromOrbital++){
											if (abs(tCallBack.getHoppingAmplitude({kx, ky, toLayer, toX, toY, toSite, toOrbital, 0}, {kx, ky, fromLayer, fromX,fromY,fromSite, fromOrbital, 0})) > EPSILON){
												for (int spin = 0; spin<2; spin++){
													if(toLayer == fromLayer && toX == fromX && toY == fromY && toSite == fromSite)
														continue;
													model << HoppingAmplitude(tCallBack, {kx, ky, toLayer, toX, toY, toSite, toOrbital, spin}, {kx, ky, fromLayer, fromX,fromY,fromSite, fromOrbital, spin});// + HC;										
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

}


//RealSpace - LatticeInformation
Geometry Graphene::setupGeometry(const vector<int> &kPoint){
	switch(grapheneUnitCell){
	
    case SelectGrapheneUnitCell::Rectangular:
		return setupRectangularGeometry(kPoint);
	
    case SelectGrapheneUnitCell::Hexagonal:
		return setupHexagonalGeometry(kPoint);
	
    default:
		TBTKExit(
			"TCallBack::radialAndAngularDependence",
			"Unknown mode.",
			"This should never happen."
		);
	}
}

Geometry Graphene::setupHexagonalGeometry(const vector<int> &kPoint){

    TBTKExit(
        "SetupHexagonalGeometry()",
        "Still building the function",
        "");
    exit(1);
}

//OK
Geometry Graphene::setupRectangularGeometry(const vector<int> &kPoint){
    Timer::tick("Setup geometry");
    vector<vector<double>> deformation = {{1, 1},{parameters.SIZE_X[0]/((double)parameters.SIZE_X[1]), parameters.SIZE_Y[0]/((double)parameters.SIZE_Y[1])}};
    
    double ap = parameters.a/sqrt(3); 
    double latticeVector[2]={3*ap, parameters.a}; 
    
    double carbonLayerSeparation = 3.3;
    Vector3d R_layer({0,0, carbonLayerSeparation});
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

    if (parameters.abStacking==true){
        Vector3d layerTranslation[2]={
            Vector3d({ 0, 0, 0}),
            Vector3d({ap, 0, 0})
        };

        for (int layer=0; layer<2; layer++){
            for (int site=0; site<parameters.numAtomsUnitCell; site++){
                offSetSite[layer][site] = offSetSite[layer][site] + layerTranslation[layer];
            }
        }
    }

    for (int layer=0; layer<2; layer++){
        for (int site=0; site<parameters.numAtomsUnitCell; site++){
            offSetSite[layer][site].x *= deformation[layer][0];
            offSetSite[layer][site].y *= deformation[layer][1];
        }
    }

    Geometry& geometry = model.getGeometry(); 
	
	int kx = kPoint[0];
	int ky = kPoint[1];
	for(int layer=0; layer<2; layer++){
		Vector3d layerPosition = R_layer;
		for(int x = 0; x < parameters.SIZE_X[layer]; x++){
			Vector3d xVector({x*latticeVector[0],0,0});
			for (int y = 0; y < parameters.SIZE_Y[layer]; y++){
				Vector3d yVector({0,y*latticeVector[1],0});
                for (int orbital=0; orbital<4; orbital++){
                    for(int spin=0; spin<2; spin++){
                        //Creates a vector to point to the current plaquette. 
                        Vector3d plaquettePosition = layerPosition+xVector*deformation[layer][0]+yVector*deformation[layer][1];	
                        geometry.setCoordinate({kx, ky, layer, x, y, 0, orbital, spin}, (plaquettePosition+offSetSite[layer][0]).getStdVector());
                        geometry.setCoordinate({kx, ky, layer, x, y, 1, orbital, spin}, (plaquettePosition+offSetSite[layer][1]).getStdVector());
                        geometry.setCoordinate({kx, ky, layer, x, y, 2, orbital, spin}, (plaquettePosition+offSetSite[layer][2]).getStdVector());
                        geometry.setCoordinate({kx, ky, layer, x, y, 3, orbital, spin}, (plaquettePosition+offSetSite[layer][3]).getStdVector());
                    }
                }                
			}
		}
	}

    Timer::tock();
    return geometry; 
}
//OK
void Graphene::printGeometry(){
    Timer::tick("Print geometry");
    ofstream fout0("systemInfo/geometry/coordinatesLayer0");
    ofstream fout1("systemInfo/geometry/coordinatesLayer1");
    ofstream fout2("systemInfo/geometry/basisVector");
    
	Geometry& geometry = model.getGeometry(); 
	for(int layer=0; layer<2; layer++){				
		for(int x = 0; x < parameters.SIZE_X[layer]; x++){
			for (int y = 0; y < parameters.SIZE_Y[layer]; y++){						
				for(int site = 0; site < parameters.numAtomsUnitCell; site++){
					Vector3d position = geometry.getCoordinate({0, 0, layer, x, y, site, 0, 0});
					if (layer == 0)	
						fout0 << position.x << "\t" << position.y << "\n";
					else
						fout1 << position.x << "\t" << position.y << "\n";
				}
			}
		}
	}

    fout0.close();
    fout1.close();
    fout2.close();
    Timer::tock();

}

//Run the Model
//OK
void Graphene::setupAndRunSolver(){
    //Setup and run Solver
    Timer::tick("Run the solver");
    // Solver::BlockDiagonalizer solver;
    solver.setModel(model);
    solver.run();
    Timer::tock();
}

//BandStructure
//TO CHECK
void Graphene::runBandStructureCalculation(){
    vector<vector<int>> kPoints = kPathsToKpoints(generateKPaths());
    parameters.numKpoints = kPoints.size();

    for (unsigned int k=0; k<parameters.numKpoints; k++){
        model = Model();
        setupGeometry(kPoints[k]);
        if (k==0){
            printGeometry();
        }
        setupModel(kPoints[k]);
        setupAndRunSolver();
        calculateBandStructure(kPoints[k], k);
    }
    plotBandStructure();    
}
//OK
void Graphene::calculateBandStructure(const std::vector<int> &kPoint, unsigned int linearKIndex){
    Timer::tick("Calculate BandStructure");

    PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

    unsigned int numBands = model.getBasisSize();

    if(bandStructure.getRanges().size()==0)//initialize the variable 'bandStructure'
   		bandStructure = Array<double>({numBands, 5*parameters.K_POINTS_PER_PATH}, 0);

    for (unsigned int band=0; band < numBands; band++){
		bandStructure[{band, linearKIndex}] = propertyExtractor.getEigenValue(band);
	}		

	Timer::tock();
}
//OK
void Graphene::plotBandStructure(){

	unsigned int numBands = model.getBasisSize();

    double min = bandStructure[{0,0}];
    double max = bandStructure[{numBands-1,0}];

    for(unsigned int n=0; n<((unsigned int) 5)*parameters.K_POINTS_PER_PATH; n++){
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
    plotter.setBoundsY(-5,5);
    for (int band=0; band<(int)numBands; band++){
        plotter.plot(bandStructure.getSlice({band, _a_}), {{"color", "black"}, {"linestyle", "-"}});
    }
    for (unsigned int n=0; n<5; n++){
        plotter.plot({((double) n)*parameters.K_POINTS_PER_PATH, ((double) n)*parameters.K_POINTS_PER_PATH},{min, max},{{"color", "black"},{"linestyle","-"}});
    }
    plotter.save("figures/bandStructure.png");
    plotter.clear();
}

//DensityOfState
//TO CHECK
void Graphene::runDOSCalculation(){
        vector<vector<int>> kPoints = generateAllKPoints();
        
        Property::DOS dos;
        for(unsigned int k=0; k<kPoints.size(); k++){
            model = Model();
            setupGeometry(kPoints[k]);
            if(k==0){
                printGeometry();
            }
            setupAndRunSolver();
            if(k==0)
                dos = calculateDOS();
            else
                dos += calculateDOS();
        }
        plotDOS(dos);
    }
//OK
Property::DOS Graphene::calculateDOS(){
    Timer::tick("Calculate DOS");

    //Property Extractor
    PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
    //Setting the energyWindow
    propertyExtractor.setEnergyWindow(parameters.LOWER_BOUND, parameters.UPPER_BOUND, parameters.RESOLUTION);
    Property::DOS dos = propertyExtractor.calculateDOS();
    Timer::tock();
    return dos;
} 
//OK
void Graphene::plotDOS(Property::DOS dos){
    Timer::tick("Plot DOS");
    dos = Smooth::gaussian(dos, parameters.SMOOTHING_SIGMA, parameters.SMOOTHING_WINDOW);
    //Plotting
    Plotter plotter;
    plotter.plot(dos);
    plotter.save("figures/DOS.png");
    plotter.clear();
    Timer::tock();
}