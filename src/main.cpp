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

int main(int argc, char **argv){
	//Initialize TBTK.
	// UnitHandler::setScales({"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Lattice parameters
	double a = 2.5;
	const int SIZE_X = 10;
	const int SIZE_Y = 10;

	// //Hamiltonian parameters;
	// complex<double> t=1.0;
	// double gamma = 0.3; 
	// double voltage = 0.0;

	// //Auxiliar
	// int SITE0 = 0;
	// int SITE1 = 1;
	// int SITE2 = 2;
	// int SITE3 = 3;

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
	Vector3d reciprocalVector[3];
	for(unsigned int n=0; n<3; n++){
		reciprocalVector[n] = 2*M_PI*latticeVector[(n+1)%3]*latticeVector[(n+2)%3]/(
			Vector3d::dotProduct(latticeVector[n], latticeVector[(n+1)%3]*latticeVector[(n+2)%3])
		);
	}

	//Setup BZ
	//Como é este vetor? Tem tamanho 10 para x, tamanho 10 para y? tipo [],[],[],[] .... []
	vector<unsigned int> numMeshPoints = {SIZE_X, SIZE_Y};
	BrillouinZone brillouinZone({
		{reciprocalVector[0].x, reciprocalVector[0].y},{reciprocalVector[1].x, reciprocalVector[1].y}
		},
		SpacePartition::MeshType::Nodal
	);
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(numMeshPoints);

	return 0;
}
