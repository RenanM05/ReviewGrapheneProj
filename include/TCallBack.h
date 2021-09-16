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
	TCallBack(double a, const Model& model): a(a), model(model){
	}

	complex<double> getHoppingAmplitude(const Index &to, const Index &from) const{	
		
		int kx = to[0];
		int ky = to[1];
			
		int toLayer = to[2];
		int fromLayer = from[2];
		
		int toX = to[3];
		int fromX = from[3];
		
		int toY = to[4];
		int fromY = from[4];

		if(toLayer == fromLayer){
			int layer = toLayer;
						
			if (toX == fromX && toY == fromY){ 
				return -t*f(to, from);
			}

			// To connect plaquette
			if (toY == fromY){
				if (toX == SIZE_X[layer]-1 && fromX == 0){
					double KX = ((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2)));
					return -t*f(to, from)*exp(i*KX*(unitCellSize[0]));
				} else if (fromX == SIZE_X[layer]-1 && toX == 0){
					double KX = ((M_PI)/(unitCellSize[0]))*(kx/((double)(SIZE_KX/2)));
					return -t*f(to, from)*exp(-i*KX*(unitCellSize[0]));
				} else{
					return -t*f(to, from);
				}
				
			} else if(toX == fromX){
				if (toY == SIZE_Y[layer]-1 && fromY == 0){
					double KY = ((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2)));
					return -t*f(to, from)*exp(i*KY*(unitCellSize[1]));
				} else if (fromY == SIZE_Y[layer]-1 && toY == 0){
					double KY = ((M_PI)/(unitCellSize[1]))*(ky/((double)(SIZE_KY/2)));
					return -t*f(to, from)*exp(-i*KY*(unitCellSize[1]));
				} else {
					return -t*f(to, from);
				}				
			} else{
				TBTKExit(
					"TCallBack::getHoppingAmplitude()",
					"Should never happened",
					""
				);
			} 
		}
		else{
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
			TBTKExit(
				"TCallBack::getHoppingAmplitude()",
				"Should never happened",
				""
			);
		} 
	}

	Vector3d getDistanceMinimizingTranslation(const Vector3d &toCoordinate, const Vector3d &fromCoordinate) const {
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
	complex<double> f(const Index &to, const Index &from) const{
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

	void setT(complex<double> t){this->t = t;}
	void setSIZE_KX(int SIZE_KX){this-> SIZE_KX = SIZE_KX;}
	void setSIZE_KY(int SIZE_KY){this-> SIZE_KY = SIZE_KY;}
	void setSIZE_X(vector<int> SIZE_X){this-> SIZE_X = SIZE_X;}
	void setSIZE_Y(vector<int> SIZE_Y){this-> SIZE_Y = SIZE_Y;}
	void setUnitCellSize(vector<double> unitCellSize){this-> unitCellSize = unitCellSize;}
	void setUnitCellBasis(vector<Vector3d> unitCellBasis){this-> unitCellBasis = unitCellBasis;}

private:
	complex<double> t;
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
};

#endif