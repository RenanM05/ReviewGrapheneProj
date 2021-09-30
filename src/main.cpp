#include "TBTK/TBTK.h"
#include "TBTK/Range.h"
#include "TBTK/Model.h"
#include "TBTK/BrillouinZone.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Smooth.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Timer.h"
#include "TBTK/TBTKMacros.h"

#include <complex>
#include <iostream>
#include <fstream>
#include <limits>

#include "TCallBack.h"
#include "Graphene.h"

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const complex<double> i(0,1);
int main(int argc, char **argv){
		//Initialize TBTK.
	Initialize();
	Timer::tick("Full Calculation");
	Timer::createAccumulator("0");
	Timer::createAccumulator("1");
	Timer::createAccumulator("2");
	Timer::createAccumulator("3");
	Timer::createAccumulator("4");
	Timer::createAccumulator("5");
	Timer::createAccumulator("6");
	Timer::createAccumulator("7");
	Timer::createAccumulator("8");
	Timer::createAccumulator("9");
	Timer::createAccumulator("10");


	//Units
	UnitHandler::setScales({"1 rad", "1 C","1 pcs","1 eV","1 Ao","1 K", "1 s"});
	Model model;

	//Initialize Graphene
	// Graphene grapheneDOS(12);
	// grapheneDOS.runDOSCalculation();

	Graphene grapheneBandStructure(36);
	grapheneBandStructure.runBandStructureCalculation();

	Timer::printAccumulators();
	Timer::tock();

	return 0;
};