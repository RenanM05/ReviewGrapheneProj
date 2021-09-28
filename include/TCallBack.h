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

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const complex<double> i(0,1);


class TCallBack: public HoppingAmplitude::AmplitudeCallback{
public: 
	enum class RadialAndAngularMode{Full, NearestNeighbor};

	TCallBack(double a, const Model& model);

	complex<double> getHoppingAmplitude(const Index &to, const Index &from) const;

	void setT(double t){this->t = t;}
	void setSIZE_KX(int SIZE_KX){this-> SIZE_KX = SIZE_KX;}
	void setSIZE_KY(int SIZE_KY){this-> SIZE_KY = SIZE_KY;}
	void setSIZE_X(vector<int> SIZE_X){this-> SIZE_X = SIZE_X;}
	void setSIZE_Y(vector<int> SIZE_Y){this-> SIZE_Y = SIZE_Y;}
	void setUnitCellSize(vector<double> unitCellSize){this-> unitCellSize = unitCellSize;}
	void setUnitCellBasis(vector<Vector3d> unitCellBasis){this-> unitCellBasis = unitCellBasis;}

private:
	double t;
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

	RadialAndAngularMode radialAndAngularMode;

	Vector3d getDistanceMinimizingTranslation(const Vector3d &toCoordinate, const Vector3d &fromCoordinate) const;
	complex<double> f(const Index &to, const Index &from) const;

	complex<double> radialAndAngularDependence(const TBTK::Index &to, const TBTK::Index &from) const;
	complex<double> radialAndAngularDependenceFull(const Index &to, const Index &from) const;
	complex<double> radialAndAngularDependenceNearestNeighbor(const Index &to, const Index &from) const;
};

#endif 