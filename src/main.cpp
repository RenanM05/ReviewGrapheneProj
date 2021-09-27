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