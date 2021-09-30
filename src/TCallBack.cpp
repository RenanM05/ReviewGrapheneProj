#include "TBTK/TBTK.h"
#include "TBTK/Range.h"
#include "TBTK/Model.h"
#include "TBTK/BrillouinZone.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Smooth.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Timer.h"
#include "TBTK/TBTKMacros.h"

#include <complex>
#include <iostream>
#include <fstream>
#include <limits>

#include "TCallBack.h"

using namespace std;
using namespace TBTK;

TCallBack::TCallBack(double a, const Model& model): model(&model), slaterKosterCalculator(a){
	this->a = a;
};

complex<double> TCallBack::getHoppingAmplitude(const Index &to, const Index &from) const{
	Timer::tick(0);
	const complex<double> i(0,1);

	int kx = to[0];
	int ky = to[1];

	const Geometry& geometry = model->getGeometry();
	const Vector3d& toCoordinate = geometry.getCoordinate(to);
	const Vector3d& fromCoordinate = geometry.getCoordinate(from);

	Vector3d k({
		((2*M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX))),
		((2*M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY))),
		0
	});
	Timer::tock(0);
	complex<double> result = 0;
	for (int x=-1; x<2; x++){
		Vector3d xTranslation = x*unitCellBasis[0];
		for (int y=-1; y<2; y++){
			Vector3d translation = xTranslation + y*unitCellBasis[1];
			result -= radialAndAngularDependence(
				to, 
				from, 
				toCoordinate, 
				fromCoordinate+translation
			)*exp(i*Vector3d::dotProduct(k, translation));
		}
	};
	return result;
} 

// This function is responsable to handle with one pair of indices
complex<double> TCallBack::radialAndAngularDependence(
	const Index &to, 
	const Index &from, 
	const Vector3d &toCoordinate,
	const Vector3d &fromCoordinate
) const{

	Vector3d difference = toCoordinate-fromCoordinate; 

	SlaterKosterCalculator::Orbital toOrbital = getOrbitalFromSubIndex(to[6]);
	SlaterKosterCalculator::Orbital fromOrbital = getOrbitalFromSubIndex(from[6]);

	return slaterKosterCalculator.calculate(toOrbital, fromOrbital, difference);
}

void TCallBack::setSIZE_KX(int SIZE_KX){
	this->SIZE_KX = SIZE_KX;
}
void TCallBack::setSIZE_KY(int SIZE_KY){
	this->SIZE_KY = SIZE_KY;
}
void TCallBack::setSIZE_X(vector<int> SIZE_X){
	this->SIZE_X = SIZE_X;
}
void TCallBack::setSIZE_Y(vector<int> SIZE_Y){
	this->SIZE_Y = SIZE_Y;
}
void TCallBack::setUnitCellSize(vector<double> unitCellSize){
	this->unitCellSize = unitCellSize;
}
void TCallBack::setUnitCellBasis(vector<Vector3d> unitCellBasis){
	this->unitCellBasis = unitCellBasis;
}