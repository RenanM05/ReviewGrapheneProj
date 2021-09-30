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

#include "Parameters.h"
#include "TCallBack.h"

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

class Graphene{
public:
	Graphene(unsigned int SIZE_K);

    void runBandStructureCalculation();
    void runDOSCalculation();

private:

	Parameters parameters;
	
	Array<double> bandStructure;

 	//Initializing the model
	Model model;
	Solver::BlockDiagonalizer solver;
 	
	//Hamiltonian Parameters.
	double t;
	TCallBack tCallBack;

	vector<vector<int>> generateAllKPoints();
	vector<vector<vector<int>>> generateKPaths();
	vector<vector<int>> kPathsToKpoints(const vector<vector<vector<int>>> &kPaths);

    void setupModel(const std::vector<int> &kPoint);
	void addHoppingAmplitude(const std::vector<int> &kPoint);
		
	TBTK::Geometry setupGeometry(const vector<int> &kPoint);
	TBTK::Geometry setupRectangularGeometry(const vector<int> &kPoint);
	TBTK::Geometry setupHexagonalGeometry(const vector<int> &kPoint);
	
	void printGeometry();
	void setupAndRunSolver();
	void calculateBandStructure(const std::vector<int> &kPoint, unsigned int linearKIndex);
	TBTK::Property::DOS calculateDOS();
	void plotDOS(TBTK::Property::DOS dos);
	void plotBandStructure();
};

#endif