#ifndef UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_GRAPHENE
#define UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_GRAPHENE

#include "TBTK/Array.h"
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

#include<TCallBack.h>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

class Graphene{
public:
	Graphene(unsigned int SIZE_K);

    void runBandStructureCalculation();
    void runDOSCalculation();

private:
	//Reciprocal-Space - BrillouinZone
	const int SIZE_K;
	const int SIZE_KX;
	const int SIZE_KY;
	unsigned int K_POINTS_PER_PATH;
	unsigned int numKpoints;
	
	Array<double> bandStructure;

 	//Initializing the model
	Model model;
	Solver::BlockDiagonalizer solver;
 	//RealSpace - LatticeInformation
	double a; 
 	const int numAtomsUnitCell;
	vector<int> SIZE_X;
	vector<int> SIZE_Y;
	bool rectangularUnitCell;
	bool addGrapheneBottomLayer;
	bool abStacking;
	bool addFluorine;
	vector<double> unitCellSize;
	vector<Vector3d> unitCellBasis;
	//Hamiltonian Parameters.
	double t;
	TCallBack tCallBack;
	//Setup energy window
	const double LOWER_BOUND;
	const double UPPER_BOUND;
	const int RESOLUTION;
	//Gaussian smoothing parameters
	const double SMOOTHING_SIGMA;
	const unsigned int SMOOTHING_WINDOW;

	vector<vector<int>> generateAllKPoints();
	vector<vector<vector<int>>> generateKPaths();
	vector<vector<int>> kPathsToKpoints(const vector<vector<vector<int>>> &kPaths);

    void setupModel(const std::vector<int> &kPoint);
	void addHoppingAmplitude(const std::vector<int> &kPoint);
	void setupGeometry(const std::vector<int> &kPoint);
	void setupRectangularGeometry(const vector<int> &kPoint);
	void setupHexagonalGeometry();
	void printGeometry();
	void setupAndRunSolver();
	void calculateBandStructure(const std::vector<int> &kPoint, unsigned int linearKIndex);
	TBTK::Property::DOS calculateDOS();
	void plotDOS(TBTK::Property::DOS dos);
	void plotBandStructure();
};

#endif