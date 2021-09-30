#ifndef UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_GRAPHENE_PARAMETERS
#define UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_GRAPHENE_PARAMETERS

#include <vector>
#include "TBTK/Vector3d.h"

class Parameters{
public:
    
    //Reciprocal-Space - BrillouinZone
	const int SIZE_K;
	const int SIZE_KX;
	const int SIZE_KY;
	unsigned int K_POINTS_PER_PATH;
	unsigned int numKpoints;
	
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
    
    //Setup energy window
	const double LOWER_BOUND;
	const double UPPER_BOUND;
	const int RESOLUTION;
	//Gaussian smoothing parameters
	const double SMOOTHING_SIGMA;
	const unsigned int SMOOTHING_WINDOW;

	enum class SelectGrapheneUnitCell{Hexagonal, Rectangular};
	SelectGrapheneUnitCell grapheneUnitCell;
	
    Parameters(unsigned int SIZE_K);

};


#endif