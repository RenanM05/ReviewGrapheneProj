#ifndef UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_SLATER_KOSTER_CALCULATOR
#define UPPSALA_RENAN_GRAPHENE_MOIRRE_PROJECT_REVIEW_SLATER_KOSTER_CALCULATOR

#include "TBTK/TBTK.h"
#include "TBTK/Timer.h"
#include "TBTK/TBTKMacros.h"
#include "TBTK/Vector3d.h"

#include <complex>

class SlaterKosterCalculator{
public:
    enum class Orbital{s,px,py,pz};
    SlaterKosterCalculator(double a);
    std::complex<double> calculate(Orbital toOrbital, Orbital fromOrbital, const TBTK::Vector3d& separation) const;

private:

    double a; 
    double t;
    double beta;
    double alpha; 

    double Vsssigma;
    double Vppsigma;
    double Vspsigma;
    double Vpppi;
    
    std::complex<double> calculateSS( const TBTK::Vector3d &separation) const;
    std::complex<double> calculateSPX(const TBTK::Vector3d &separation) const;
    std::complex<double> calculateSPY(const TBTK::Vector3d &separation) const;
    std::complex<double> calculateSPZ(const TBTK::Vector3d &separation) const;

    std::complex<double> calculatePXS( const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePXPX(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePXPY(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePXPZ(const TBTK::Vector3d &separation) const;

    std::complex<double> calculatePYS(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePYPX(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePYPY(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePYPZ(const TBTK::Vector3d &separation) const;

    std::complex<double> calculatePZS(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePZPX(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePZPY(const TBTK::Vector3d &separation) const;
    std::complex<double> calculatePZPZ(const TBTK::Vector3d &separation) const;
};

#endif