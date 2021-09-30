#include "TBTK/Model.h"
#include "TBTK/TBTK.h"
#include "TBTK/Model.h"

#include <complex>
#include <iostream>
#include <fstream>

#include "SetGeometry.h"
#include "Graphene.h"

using namespace std;
using namespace TBTK;

SetGeometry::SetGeometry():
//RealSpace - LatticeInformation
a(2.45),
numAtomsUnitCell(4),
SIZE_X({2,2}),
SIZE_Y({2,2}),
rectangularUnitCell(true),
addGrapheneBottomLayer(false),
abStacking(true),
addFluorine(false),
unitCellSize({SIZE_X[0]*a*sqrt(3), SIZE_Y[0]*a}),
unitCellBasis({{unitCellSize[0], 0, 0},{0, unitCellSize[1], 0}})
{
    grapheneUnitCell = SelectGrapheneUnitCell::Rectangular;
}


Geometry SetGeometry::setupGeometry(const vector<int> &kPoint){
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

Geometry SetGeometry::setupHexagonalGeometry(const vector<int> &kPoint){

    TBTKExit(
        "SetupHexagonalGeometry()",
        "Still building the function",
        "");
    exit(1);
}


Geometry SetGeometry::setupRectangularGeometry(const vector<int> &kPoint){
    Timer::tick("Setup geometry");
    vector<vector<double>> deformation = {{1, 1},{SIZE_X[0]/((double)SIZE_X[1]), SIZE_Y[0]/((double)SIZE_Y[1])}};
    
    double ap = a/sqrt(3); 
    double latticeVector[2]={3*ap, a}; 
    
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
	
	int kx = kPoint[0];
	int ky = kPoint[1];
	for(int layer=0; layer<2; layer++){
		Vector3d layerPosition = R_layer;
		for(int x = 0; x < SIZE_X[layer]; x++){
			Vector3d xVector({x*latticeVector[0],0,0});
			for (int y = 0; y < SIZE_Y[layer]; y++){
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
void SetGeometry::printGeometry(){
    Timer::tick("Print geometry");
    ofstream fout0("systemInfo/geometry/coordinatesLayer0");
    ofstream fout1("systemInfo/geometry/coordinatesLayer1");
    ofstream fout2("systemInfo/geometry/basisVector");
    
	Geometry& geometry = model.getGeometry(); 
	for(int layer=0; layer<2; layer++){				
		for(int x = 0; x < SIZE_X[layer]; x++){
			for (int y = 0; y < SIZE_Y[layer]; y++){						
				for(int site = 0; site < numAtomsUnitCell; site++){
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
