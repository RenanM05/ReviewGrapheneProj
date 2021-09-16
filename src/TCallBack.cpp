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

TCallBack::TCallBack(double a, const Model& model): a(a), model(model){
}

complex<double> TCallBack::getHoppingAmplitude(const Index &to, const Index &from) const{	
    
    int kx = to[0];
    int ky = to[1];
    
    const Geometry& geometry = model.getGeometry();
    Vector3d toCoordinate = geometry.getCoordinate(to);
    Vector3d fromCoordinate = geometry.getCoordinate(from);
    Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);

    Vector3d k({
        ((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2))),
        ((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2))),
        0
    });
    
    return -t*f(to, from)*exp(i*Vector3d::dotProduct(k, distanceMinimizingTranslation));
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

// This function is responsable to handle with one pair of indices
complex<double> TCallBack::f(const Index &to, const Index &from) const{
    Vector3d toCoordinate=model.getGeometry().getCoordinate(to);
    Vector3d fromCoordinate=model.getGeometry().getCoordinate(from);

    Vector3d distanceMinimizingTranslation = getDistanceMinimizingTranslation(toCoordinate, fromCoordinate);
    Vector3d difference = fromCoordinate + distanceMinimizingTranslation - toCoordinate;
    double smallestDistance = difference.norm();

    double beta=3; 		

    if (to[2] == from[2]){
        return exp(-beta*(smallestDistance/(a/sqrt(3)) - 1));
    }
    else{
        return Vector3d::dotProduct(difference.unit(), {0, 0, 1})*exp(-beta*(smallestDistance/(a/sqrt(3)) - 1));
    }
}

void TCallBack::setT(complex<double> t){this->t = t;}
void TCallBack::setSIZE_KX(int SIZE_KX){this-> SIZE_KX = SIZE_KX;}
void TCallBack::setSIZE_KY(int SIZE_KY){this-> SIZE_KY = SIZE_KY;}
void TCallBack::setSIZE_X(vector<int> SIZE_X){this-> SIZE_X = SIZE_X;}
void TCallBack::setSIZE_Y(vector<int> SIZE_Y){this-> SIZE_Y = SIZE_Y;}
void TCallBack::setUnitCellSize(vector<double> unitCellSize){this-> unitCellSize = unitCellSize;}
void TCallBack::setUnitCellBasis(vector<Vector3d> unitCellBasis){this-> unitCellBasis = unitCellBasis;}