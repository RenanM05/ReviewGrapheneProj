#include "TBTK/TBTK.h"
#include "SlaterKosterCalculator.h"

using namespace std;
using namespace TBTK;

SlaterKosterCalculator::SlaterKosterCalculator(double a):
a(a),
t(2.8),
beta(3),//1.9467973982212834;
alpha(-0.2*exp(+beta*(3.4/(a/sqrt(3)) - 1))),
// alpha(-0.4*exp(+beta*(3.4/(a/sqrt(3)) - 1))),
Vsssigma(1),
Vppsigma(alpha),
Vspsigma(1),
Vpppi(-t)
{
}

complex<double> SlaterKosterCalculator::calculate(Orbital toOrbital, Orbital fromOrbital, const TBTK::Vector3d& separation) const{
    if (separation.norm()<1e-3)
        return 0;
 
    if(toOrbital==Orbital::s && fromOrbital==Orbital::s){
        return calculateSS(separation);
    } else if (toOrbital==Orbital::s && fromOrbital==Orbital::px){
        return calculateSPX(separation);
    } else if (toOrbital==Orbital::s && fromOrbital==Orbital::py){
        return calculateSPY(separation); 
    }else if  (toOrbital==Orbital::s && fromOrbital==Orbital::pz){
        return calculateSPY(separation);  
    }

    if(toOrbital==Orbital::px && fromOrbital==Orbital::px){
        return calculatePXPX(separation);
    } else if (toOrbital==Orbital::px && fromOrbital==Orbital::s){
        return calculatePXS(separation);
    } else if (toOrbital==Orbital::px && fromOrbital==Orbital::py){
        return calculatePXPY(separation); 
    }else if  (toOrbital==Orbital::px && fromOrbital==Orbital::pz){
        return calculatePXPZ(separation);  
    }

    if(toOrbital==Orbital::py && fromOrbital==Orbital::py){
        return calculatePYPY(separation);
    } else if (toOrbital==Orbital::py && fromOrbital==Orbital::s){
        return calculatePYS(separation);
    } else if (toOrbital==Orbital::py && fromOrbital==Orbital::px){
        return calculatePYPX(separation); 
    }else if  (toOrbital==Orbital::py && fromOrbital==Orbital::pz){
        return calculatePYPZ(separation);  
    }

    if(toOrbital==Orbital::pz && fromOrbital==Orbital::pz){
        return calculatePZPZ(separation);
    } else if (toOrbital==Orbital::pz && fromOrbital==Orbital::s){
        return calculatePZS(separation);
    } else if (toOrbital==Orbital::pz && fromOrbital==Orbital::px){
        return calculatePZPX(separation); 
    }else if  (toOrbital==Orbital::pz && fromOrbital==Orbital::py){
        return calculatePZPY(separation);  
    }

    TBTKExit(
    "SlaterKosterCalculator::calculate()",
    "Unknown orbital",
    "This should never happen."
    );
}

//S
complex<double> SlaterKosterCalculator::calculateSS(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    return Vsssigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculateSPX(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    return l*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculateSPY(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return m*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculateSPZ(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    return n*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}

//PX
complex<double> SlaterKosterCalculator::calculatePXS(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    return l*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePXPX(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    return ((l*l)*Vppsigma+(1-(l*l))*Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePXPY(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return l*m*(Vppsigma - Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePXPZ(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    return l*n*(Vppsigma - Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}


//PY
complex<double> SlaterKosterCalculator::calculatePYS(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return m*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePYPX(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return l*m*(Vppsigma - Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePYPY(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return ((m*m)*Vppsigma+(1-(m*m))*Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePYPZ(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    return m*n*(Vppsigma - Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
//PZ
complex<double> SlaterKosterCalculator::calculatePZS(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    return n*Vspsigma*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePZPX(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
    double l = Vector3d::dotProduct(separation.unit(), {1, 0, 0});
    double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    return l*n*(Vppsigma-Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePZPY(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
	double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
    double m = Vector3d::dotProduct(separation.unit(), {0, 1, 0});
    return n*m*(Vppsigma-Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
}
complex<double> SlaterKosterCalculator::calculatePZPZ(const TBTK::Vector3d& separation) const{
    double distance = separation.norm();
	double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});
	return ((n*n)*Vppsigma+(1-(n*n))*Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
} 

//     double distance = separation.norm();
// 	double n = Vector3d::dotProduct(separation.unit(), {0, 0, 1});

// 	double beta=1.9467973982212834;
// 	double alpha=-0.2*exp(+beta*(3.3/(a/sqrt(3)) - 1));
// 	double Vppsigma = alpha*exp(-beta*(distance/(a/sqrt(3)) - 1));
// 	double Vpppi = -t*exp(-beta*(distance/(a/sqrt(3)) - 1));

// 	return (n*n)*Vppsigma+(1-(n*n))*Vpppi;

// 	// return ((n*n)*Vppsigma+(1-(n*n))*Vpppi)*exp(-beta*(distance/(a/sqrt(3)) - 1));
// }