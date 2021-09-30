#ifndef UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_T_CALL_BACK
#define UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_T_CALL_BACK

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

#include "SlaterKosterCalculator.h"

class TCallBack: public TBTK::HoppingAmplitude::AmplitudeCallback{
public: 
	TCallBack(double a, const TBTK::Model& model);
	std::complex<double> getHoppingAmplitude(const TBTK::Index &to, const TBTK::Index &from) const;

	void setSIZE_KX(int SIZE_KX);
	void setSIZE_KY(int SIZE_KY);
	void setSIZE_X(std::vector<int> SIZE_X);
	void setSIZE_Y(std::vector<int> SIZE_Y);
	void setUnitCellSize(std::vector<double> unitCellSize);
	void setUnitCellBasis(std::vector<TBTK::Vector3d> unitCellBasis);

private:

	std::vector<int> SIZE_X;
	std::vector<int> SIZE_Y;
	std::vector<double> unitCellSize;
	std::vector<TBTK::Vector3d> unitCellBasis;

	int SIZE_KX;
	int SIZE_KY;
	double R_X;
	double R_Y;
	double a;
	const TBTK::Model *model;

	SlaterKosterCalculator slaterKosterCalculator;
	
	std::complex<double> radialAndAngularDependence(
		const TBTK::Index &to,
		const TBTK::Index &from,
		const TBTK::Vector3d &toCoordinate,
		const TBTK::Vector3d &fromCoordinate
	) const; 
	
	SlaterKosterCalculator::Orbital getOrbitalFromSubIndex(int subIndex) const;
};

inline SlaterKosterCalculator::Orbital TCallBack::getOrbitalFromSubIndex(int subIndex) const{
	switch(subIndex){
	case 0:
		return SlaterKosterCalculator::Orbital::s;
	case 1:
		return SlaterKosterCalculator::Orbital::px;
	case 2:
		return SlaterKosterCalculator::Orbital::py;
	case 3:
		return SlaterKosterCalculator::Orbital::pz;
	default:
		TBTKExit(
			"TCallBack::getOrbitalSubIndex()",
			"Invalid subindex '" << subIndex << "'. Must be 0-3.",
			""
		);
	}
}

#endif 