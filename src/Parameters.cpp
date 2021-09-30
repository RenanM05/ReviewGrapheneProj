#include "Parameters.h"

Parameters::Parameters(unsigned int SIZE_K):
    //Reciprocal-Space - BrillouinZone
    SIZE_K(SIZE_K),
    SIZE_KX(SIZE_K),
    SIZE_KY(SIZE_K),
    K_POINTS_PER_PATH(SIZE_K/2),
    //RealSpace - LatticeInformation
    a(2.45),
    numAtomsUnitCell(4),
    SIZE_X({2,2}),
    SIZE_Y({2,2}),
    rectangularUnitCell(true),
    addGrapheneBottomLayer(false),
    abStacking(true),
    addFluorine(false),
    unitCellSize({SIZE_X[0]*a*sqrt(3), SIZE_Y[0]*a}),
    unitCellBasis({{unitCellSize[0], 0, 0},{0, unitCellSize[1], 0}}),
    //Setup energy window
    LOWER_BOUND(-2),
    UPPER_BOUND(+2),
    RESOLUTION(1000),
    //Gaussian smoothing parameters
    SMOOTHING_SIGMA(0.1),
    SMOOTHING_WINDOW(51)
{
    grapheneUnitCell = SelectGrapheneUnitCell::Rectangular;
}