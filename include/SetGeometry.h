#ifndef UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_SET_GEOMETRY
#define UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_SET_GEOMETRY

#include "TBTK/TBTK.h"
#include "Graphene.h"

class SetGeometry{
public:
    
	//RealSpace - LatticeInformation
	double a; 
 	const int numAtomsUnitCell;
	std::vector<int> SIZE_X;
	std::vector<int> SIZE_Y;
	bool rectangularUnitCell;
	bool addGrapheneBottomLayer;
	bool abStacking;
	bool addFluorine;
	std::vector<double> unitCellSize;
	std::vector<TBTK::Vector3d> unitCellBasis;

	TBTK::Geometry setupGeometry(const std::vector<int> &kPoint);
	TBTK::Geometry setupRectangularGeometry(const std::vector<int> &kPoint);
	TBTK::Geometry setupHexagonalGeometry(const std::vector<int> &kPoint);
    void printGeometry();

	SetGeometry();

private:
	//Initializing the model
	TBTK::Model model;
	
    enum class SelectGrapheneUnitCell{Hexagonal, Rectangular};
	SelectGrapheneUnitCell grapheneUnitCell;

};

#endif