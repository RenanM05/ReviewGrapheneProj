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

TCallBack::TCallBack(double a, const Model& model): a(a), model(&model)
{
    radialAndAngularMode = RadialAndAngularMode::Full;    
}

complex<double> TCallBack::getHoppingAmplitude(const Index &to, const Index &from) const{
	const complex<double> i(0,1);
	int kx = to[0];
	int ky = to[1];
	const Geometry& geometry = model->getGeometry();
	Vector3d toCoordinate = geometry.getCoordinate(to);
	Vector3d fromCoordinate = geometry.getCoordinate(from);
	Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);
	Vector3d k({
		(2*(M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX))),
		(2*(M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY))),
		0
	});
	return radialAndAngularDependence(to, from)*exp(i*Vector3d::dotProduct(k, distanceMinimizingTranslation));
} 
Vector3d TCallBack::getDistanceMinimizingTranslation(const Vector3d &toCoordinate, const Vector3d &fromCoordinate) const {
	double smallestNorm = std::numeric_limits<double>::infinity(); 		
	Vector3d requiredTranslation({0,0,0});
	for (int x=-1; x<2; x++){
		for (int y=-1; y<2; y++){
			Vector3d translation = x*unitCellBasis[0]+ y*unitCellBasis[1];
			Vector3d translatedFromCoordinate = fromCoordinate + translation;
			Vector3d difference = toCoordinate - translatedFromCoordinate;
			double norm = difference.norm();
			if(smallestNorm > norm){
				smallestNorm = norm;
				requiredTranslation = translation;
			}
		};
	};
	return requiredTranslation;
}

complex<double> TCallBack::radialAndAngularDependence(const Index &to, const Index &from) const{
	switch(radialAndAngularMode){
	case RadialAndAngularMode::Full:
		return radialAndAngularDependenceFull(to, from);
	case RadialAndAngularMode::NearestNeighbor:
		return radialAndAngularDependenceNearestNeighbor(to, from);
	default:
		TBTKExit(
			"TCallBack::radialAndAngularDependence",
			"Unknown mode.",
			"This should never happen."
		);
	}
}
complex<double> TCallBack::radialAndAngularDependenceFull(const Index &to, const Index &from) const{
	
	Vector3d toCoordinate=model->getGeometry().getCoordinate(to);
	Vector3d fromCoordinate=model->getGeometry().getCoordinate(from);

	Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);
	Vector3d difference = fromCoordinate + distanceMinimizingTranslation - toCoordinate;
	double distance = difference.norm();

	// beta: calculated analitically to give (0.2/t, t=2.8eV) for smallestdistance = 3.4
	double beta=1.9467973982212834;

	double alpha=-0.2*exp(+beta*(3.3/(a/sqrt(3)) - 1));
	double n = Vector3d::dotProduct(difference.unit(), {0, 0, 1});
	double Vppsigma = alpha*exp(-beta*(distance/(a/sqrt(3)) - 1));
	double Vpppi = -t*exp(-beta*(distance/(a/sqrt(3)) - 1));

	//Just for Pz-Pz
	return (n*n)*Vppsigma+(1-(n*n))*Vpppi;
}
complex<double> TCallBack::radialAndAngularDependenceNearestNeighbor(const Index &to, const Index &from) const{
	Vector3d toCoordinate=model->getGeometry().getCoordinate(to);
	Vector3d fromCoordinate=model->getGeometry().getCoordinate(from);

	Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);
	Vector3d difference = fromCoordinate + distanceMinimizingTranslation - toCoordinate;
	double smallestDistance = difference.norm();
	double beta=200;
    if (to[2] == from[2]){
		return -t*exp(-beta*(smallestDistance/(a/sqrt(3)) - 1));
	} else if((to[2] == 0 && from[2] == 1) || (to[2] == 1 && from[2] == 0)){
		return -t*(0.2/2.8)*exp(-beta*(smallestDistance/(3.4) - 1));
	} else if((to[2] == 1 && from[2] == 2) || (to[2] == 2 && from[2] == 1)){
		return -t*(0.2/2.8)*exp(-beta*(smallestDistance/(3.4) - 1));
	}else{
		TBTKExit(
			"TCallBack::radialAndAngularDependenceNearestNeighbor()",
			"Unexpected layer combination.",
			"This should never happen."
		);
	}	
}
