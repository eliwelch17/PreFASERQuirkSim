#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <iomanip>
#include <string>
#include <chrono>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>
#include <filesystem> 
#include <limits>
//source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc11-opt/setup.sh if on lxplus
//g++ -o quirk_run quirk_run.cxx
// Constants (based on the Mathematica script)
const double ZCu = 29;
const double ZoACu = ZCu / 63.546;
const double rhoCu = 8.960;
const double I0Cu = 322.0e-6;
const double aCu = 0.14339;
const double kCu = 2.9044;
const double x0Cu = -0.0254;
const double x1Cu = 3.2792;
const double CbarCu = 4.4190;
const double d0Cu = 0.08;

const double ZCc = 8.56;
const double ZoACc = ZCc / 40.04;
const double rhoCc = 2.300;
const double I0Cc = 135.2e-6;
const double aCc = 0.07515;
const double kCc = 3.5467;
const double x0Cc = 0.1301;
const double x1Cc = 3.0466;
const double CbarCc = 3.9464;
const double d0Cc = 0.00;

const double ZRock = 11;
const double ZoARock = ZRock / 22.99;
const double rhoRock = 2.650;
const double I0Rock = 136.4e-6;
const double aRock = 0.08301;
const double kRock = 3.4120;
const double x0Rock = 0.0492;
const double x1Rock = 3.0549;
const double CbarRock = 3.7738;
const double d0Rock = 0.00;


using namespace std;

std::vector<double> Cross(const std::vector<double>& v1, const std::vector<double>& v2) {

    //cross product function
    if (v1.size() != 3 || v2.size() != 3) {
        throw std::invalid_argument("Both vectors must be 3-dimensional");
    }
    std::vector<double> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}



std::vector<long double> CrossLong(const std::vector<long double>& v1, const std::vector<long double>& v2) {
    //cross product function with long's
    if (v1.size() != 3 || v2.size() != 3) {
        throw std::invalid_argument("Both vectors must be 3-dimensional");
    }
    std::vector<long double> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}


double Gamma(double beta) {
    return 1.0 / sqrt(1.0 - beta * beta);
}

double Emax(double m, double beta) {
    return (2 * 0.511 * beta * beta * pow(Gamma(beta), 2)) / 
           (1 + 2 * Gamma(beta) * (0.511e-3 / m) + pow(0.511e-3 / m, 2));
}

double Xi(double z, double ZoA, double rho, double beta) {
    return 0.5 * 0.307075 * z * z * ZoA * rho / (beta * beta);
}

double Delta(double x, double x0, double x1, double Cbar, double a, double k, double d0) {
    if (x >= x1) {
        return 2 * log(10) * x - Cbar;
    } else if (x >= x0) {
        return 2 * log(10) * x - Cbar + a * pow(x1 - x, k);
    } else {
        return d0 * pow(10, 2 * (x - x0));
    }
}

double Eox(double m, double z, double ZoA, double rho, double I0, double x0, double x1, double Cbar, double a, double k, double d0, double beta) {
   
    return 2 * Xi(z, ZoA, rho, beta) * (0.5 * log((2 * 0.511 * beta * beta * pow(Gamma(beta), 2) * Emax(m, beta)) / pow(I0, 2))
           - beta * beta - Delta(log10(beta * Gamma(beta)), x0, x1, Cbar, a, k, d0) / 2);
         
}


/*double EoxCuAll(double mq, int param, double beta) {
   

    // Determine the scalar value based on the magnitude
    double scalar;
    double ZCu_eff=1;
    if (beta <= 0.00226) {
        scalar = 37597.30061169589 * beta;
    } else if (beta >= 0.06345454545454546) {
        scalar = Eox(mq, ZCu_eff, ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, beta) * 0.197;
      
    } else {
        scalar = 0.34289008715835223 + 33909.5576744394 * beta + 2.6438521703577857e6 * std::pow(beta, 2)
                 - 5.582148886272709e8 * std::pow(beta, 3) + 3.916027090833348e10 * std::pow(beta, 4)
                 - 1.6320734910018296e12 * std::pow(beta, 5) + 4.608165819030902e13 * std::pow(beta, 6)
                 - 9.294312877895761e14 * std::pow(beta, 7) + 1.370510669479835e16 * std::pow(beta, 8)
                 - 1.486143034696763e17 * std::pow(beta, 9) + 1.174756123051605e18 * std::pow(beta, 10)
                 - 6.596425471576249e18 * std::pow(beta, 11) + 2.496247979943447e19 * std::pow(beta, 12)
                 - 5.714216718459019e19 * std::pow(beta, 13) + 5.982648430710585e19 * std::pow(beta, 14);
    }

    // Return a vector where each component is the same scalar value
    
    return scalar;
}*/

double EoxCuAll(double mq, int param, double beta) {
    //de/dx copper
    //EoxCuAll with recasting variables to long double
    long double mq_ld = static_cast<long double>(mq);
    long double beta_ld = static_cast<long double>(beta);
    long double scalar;
    long double ZCu_eff = 1.0L;

    if (beta_ld <= 0.00226L) {
        scalar = 37597.30061169589L * beta_ld;
    } else if (beta_ld >= 0.06345454545454546L) {
        scalar = static_cast<long double>(Eox(static_cast<double>(mq_ld), static_cast<double>(ZCu_eff), ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, static_cast<double>(beta_ld))) * 0.197L;
    } else {
        scalar = 0.34289008715835223L + 33909.5576744394L * beta_ld + 2.6438521703577857e6L * std::pow(beta_ld, 2)
                 - 5.582148886272709e8L * std::pow(beta_ld, 3) + 3.916027090833348e10L * std::pow(beta_ld, 4)
                 - 1.6320734910018296e12L * std::pow(beta_ld, 5) + 4.608165819030902e13L * std::pow(beta_ld, 6)
                 - 9.294312877895761e14L * std::pow(beta_ld, 7) + 1.370510669479835e16L * std::pow(beta_ld, 8)
                 - 1.486143034696763e17L * std::pow(beta_ld, 9) + 1.174756123051605e18L * std::pow(beta_ld, 10)
                 - 6.596425471576249e18L * std::pow(beta_ld, 11) + 2.496247979943447e19L * std::pow(beta_ld, 12)
                 - 5.714216718459019e19L * std::pow(beta_ld, 13) + 5.982648430710585e19L * std::pow(beta_ld, 14);
    }

    return static_cast<double>(scalar);
}


double EoxCcAll(double mq, int param, double beta) {
 //de/dx concrete
    double scalar;
    double ZCc_eff=1;
    if (beta <= 0.00226) {
        scalar = 30376.62714732753 * beta;
    } else if (beta >= 0.06181818181818182) {
        scalar = Eox(mq, ZCc_eff=1, ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, beta) * 0.197;
    } else {
        scalar = 0.3717416673986129 + 26377.245398907748 * beta + 2.8807310255188923e6 * std::pow(beta, 2)
                 - 6.152256215129855e8 * std::pow(beta, 3) + 4.441507306806682e10 * std::pow(beta, 4)
                 - 1.8750617184047375e12 * std::pow(beta, 5) + 5.310954532975338e13 * std::pow(beta, 6)
                 - 1.0693146651849746e15 * std::pow(beta, 7) + 1.5707786950520584e16 * std::pow(beta, 8)
                 - 1.695911162356092e17 * std::pow(beta, 9) + 1.3350858443896056e18 * std::pow(beta, 10)
                 - 7.470735431733813e18 * std::pow(beta, 11) + 2.819535817843903e19 * std::pow(beta, 12)
                 - 6.442273849826215e19 * std::pow(beta, 13) + 6.737738168445225e19 * std::pow(beta, 14);
    }

    return scalar;
}

double EoxRockAll(double mq, int param, double beta) {
   //de/dx rock
    double scalar;
    double ZRock_eff=1;
    if (beta <= 0.00226) {
        scalar = 28340.291807946152 * beta;
    } else if (beta >= 0.05745454545454545) {
        scalar = Eox(mq, ZRock_eff, ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, beta) * 0.197;
    } else {
        scalar = 0.3044358906484703 + 25065.834574193614 * beta + 2.350400267498349e6 * std::pow(beta, 2)
                 - 4.976688103007817e8 * std::pow(beta, 3) + 3.513441645396524e10 * std::pow(beta, 4)
                 - 1.462059229558537e12 * std::pow(beta, 5) + 4.113520947417071e13 * std::pow(beta, 6)
                 - 8.278484331611662e14 * std::pow(beta, 7) + 1.2214854741842604e16 * std::pow(beta, 8)
                 - 1.3296907090779789e17 * std::pow(beta, 9) + 1.0585274314396489e18 * std::pow(beta, 10)
                 - 6.003163825933505e18 * std::pow(beta, 11) + 2.30021301170228e19 * std::pow(beta, 12)
                 - 5.34284980622759e19 * std::pow(beta, 13) + 5.686153048402872e19 * std::pow(beta, 14);
    }
    return scalar;
}


//really the de/dx Gaus functions should use truncated distributions

/*double truncated_normal(double mean, double std_dev) {
    boost::math::normal_distribution<> normal_dist(mean, std_dev);
    boost::random::uniform_real_distribution<> uniform_dist(0.0, 1.0);

    double lower_cdf = boost::math::cdf(normal_dist, 0);
    double upper_cdf = boost::math::cdf(normal_dist, std::numeric_limits<double>::infinity());

    double u = uniform_dist(gen) * (upper_cdf - lower_cdf) + lower_cdf;

    // inverse CDF corresponding to the generated uniform random number
    return boost::math::quantile(normal_dist, u);
}*/



//FIXME: once validated and random is turned back on, need to pass gen as input variable!
double EoxGaus(double m, double z, double ZoA, double rho, double I0, double x0, double x1, double Cbar, double a, double k, double d0, double beta, double delta_x, std::function<double(double, int, double)> EoxAllFunc) {
    // de/dx gaus 
    double z_eff =1;
    double mean = EoxAllFunc(m, z_eff, beta);
    double std_dev = 0.197 * sqrt(Xi(z_eff, ZoA, rho, beta) * delta_x * Emax(m, beta) * (1 - beta * beta / 2)) / delta_x;
    std::normal_distribution<> d(mean, std_dev); 
    return 1.0; //FIXME: temporary

    //return d(gen); 
    //return truncated_normal(mean, std_dev); //for truncated distributions
}


// de/dx gausfor different materials
double EoxCuGaus(double mq, int param, const double beta, double dx) {
    return EoxGaus(mq, ZCu, ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, beta, dx, EoxCuAll);
}

double EoxCcGaus(double mq, int param, const double beta, double dx) {
    return EoxGaus(mq, ZCc, ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, beta, dx, EoxCcAll);
}

double EoxRockGaus(double mq, int param, const double beta, double dx) {
    return EoxGaus(mq, ZRock, ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, beta, dx, EoxRockAll);
}

std::vector<double> EoxCu(double mq, int param, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = Eox(mq, ZCu, ZoACu, rhoCu, I0Cu, x0Cu, x1Cu, CbarCu, aCu, kCu, d0Cu, v[i]);
    }
    return result;
}

std::vector<double> EoxCc(double mq, int param, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = Eox(mq, ZCc, ZoACc, rhoCc, I0Cc, x0Cc, x1Cc, CbarCc, aCc, kCc, d0Cc, v[i]);
    }
    return result;
}

std::vector<double> EoxRock(double mq, int param, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = Eox(mq, ZRock, ZoARock, rhoRock, I0Rock, x0Rock, x1Rock, CbarRock, aRock, kRock, d0Rock, v[i]);
    }
    return result;
}

int Loct(double x, double y, double z) {
    //determine location region of quirks
    if ((((sqrt(x * x + y * y) > 0.017e6) && (abs(z - 19.9e6) < 0.9e6)) ||
        ((abs(x) < (0.094 / 2) * 1e6) && (abs(y + (0.605 / 2 - 0.067) * 1e6) < (0.605 / 2) * 1e6) && (abs(z - 140.5e6) < 0.5e6)) ||
        (abs(z - 385.0e6) < 5.0e6) || (abs(z - 435.0e6) < 45.0e6))) {
        if ((sqrt(x * x + y * y) > 0.017e6) && (abs(z - 19.9e6) < 0.9e6)) return 1;
        if ((abs(x) < (0.094 / 2) * 1e6) && (abs(y + (0.605 / 2 - 0.067) * 1e6) < (0.605 / 2) * 1e6) && (abs(z - 140.5e6) < 0.5e6)) return 2;
        if (abs(z - 385.0e6) < 5.0e6) return 3;
        if (abs(z - 435.0e6) < 45.0e6) return 4;
    }
    return 0;
}

std::vector<double> Bct(double x, double y, double z) {
    //B field function of location region
    if (((sqrt(x * x + y * y) < 0.06e6) && (abs(z - 72.287e6) < 12.365e6)) ||
        ((((sqrt((x - 0.093e6) * (x - 0.093e6) + y * y) < 0.04e6) || (sqrt((x + 0.093e6) * (x + 0.093e6) + y * y) < 0.04e6)) && (abs(z - 158.2e6) < 4.725e6)))) {
        if ((sqrt(x * x + y * y) < 0.06e6) && (abs(z - 72.287e6) < 12.365e6)) return {0, 3.5, 0};
        if (((sqrt((x - 0.093e6) * (x - 0.093e6) + y * y) < 0.04e6) || (sqrt((x + 0.093e6) * (x + 0.093e6) + y * y) < 0.04e6)) && (abs(z - 158.2e6) < 4.725e6)) return {0, -3.5, 0};
    }
    return {0, 0, 0};
}

std::vector<long double> BctLong(long double x, long double y, long double z) {
    //B long function
    if (((std::sqrt(x * x + y * y) < 0.06e6L) && (std::abs(z - 72.287e6L) < 12.365e6L)) ||
        ((((std::sqrt((x - 0.093e6L) * (x - 0.093e6L) + y * y) < 0.04e6L) || 
           (std::sqrt((x + 0.093e6L) * (x + 0.093e6L) + y * y) < 0.04e6L)) && 
           (std::abs(z - 158.2e6L) < 4.725e6L)))) 
    {
        if ((std::sqrt(x * x + y * y) < 0.06e6L) && (std::abs(z - 72.287e6L) < 12.365e6L)) {
            return {0.0L, 3.5L, 0.0L};
        }
        if (((std::sqrt((x - 0.093e6L) * (x - 0.093e6L) + y * y) < 0.04e6L) || 
             (std::sqrt((x + 0.093e6L) * (x + 0.093e6L) + y * y) < 0.04e6L)) && 
             (std::abs(z - 158.2e6L) < 4.725e6L)) 
        {
            return {0.0L, -3.5L, 0.0L};
        }
    }
    return {0.0L, 0.0L, 0.0L};
}



int Layer(double x, double y, double z) {
    //determine scintilaltor layers
    if (abs(x) < 0.15e6 && abs(y) < 0.15e6 && abs(z - 480.01e6) < 0.01e6) return 1;
    if (abs(x) < 0.2e6 && abs(y) < 0.2e6 && abs(z - 481.555e6) < 0.005e6) return 2;
    if (abs(x) < 0.15e6 && abs(y) < 0.15e6 && abs(z - 484.18e6) < 0.01e6) return 3;
    return 0;
}

// Helper functions for vector operations:
std::vector<double> SubtractVectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

//std::vector<double> AddVectors(const std::vector<double>& v1, const std::vector<double>& v2) {
 //   return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
//}

std::vector<double> AddVectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result(3);
    std::vector<double> c(3, 0.0); 
    for (int i = 0; i < 3; ++i) {
        double y = v2[i] - c[i];
        double t = v1[i] + y;
        c[i] = (t - v1[i]) - y;
        result[i] = t;
    }
    return result;
}

std::vector<double> MultiplyVector(const std::vector<double>& v, double scalar) {
    return {v[0] * scalar, v[1] * scalar, v[2] * scalar};
}

std::vector<long double> MultiplyVectorLong(const std::vector<long double>& v, long double scalar) {
    return {v[0] * scalar, v[1] * scalar, v[2] * scalar};
}

//std::vector<double> DivideVector(const std::vector<double>& v, double scalar) {
 //   return {v[0] / scalar, v[1] / scalar, v[2] / scalar};
//}

std::vector<double> DivideVector(const std::vector<double>& v, double scalar) {
    if (std::abs(scalar) < std::numeric_limits<double>::epsilon()) {
        scalar = (scalar < 0) ? -std::numeric_limits<double>::epsilon() : std::numeric_limits<double>::epsilon();
    }
    
    return {v[0] / scalar, v[1] / scalar, v[2] / scalar};
}

//double DotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
//    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
//}

double DotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
    //employ Kahan Summation method
    double sum = 0.0; 
    double c = 0.0; 
    for (size_t i = 0; i < v1.size(); ++i) {
        double y = v1[i] * v2[i] - c; 
        double t = sum + y; 
        c = (t - sum) - y; 
        sum = t; 
    }
    return sum;
}

std::vector<double> Normalize(const std::vector<double>& v) {
    double norm = std::sqrt(DotProduct(v, v));
    return DivideVector(v, norm);
}

std::vector<long double> NormalizeLong(const std::vector<double>& v) {
    long double norm = std::sqrt(static_cast<long double>(DotProduct(v, v)));
    std::vector<long double> normalizedVector(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        normalizedVector[i] = static_cast<long double>(v[i]) / norm;
    }
    return normalizedVector;
}

// Function to calculate the value of ct
double CalculateCt(const std::vector<double>& v, const std::vector<double>& Beta, const std::vector<double>& r1, const std::vector<double>& r2, const std::vector<double>& F, double E1, double E2) {
    return 1 - DotProduct(v, Beta) - DotProduct(SubtractVectors(r1, r2), SubtractVectors(F, MultiplyVector(Beta, DotProduct(v, F)))) / (300000 * (E1 + E2));
}

// Function to calculate the travel distance
double CalculateDistance(const std::vector<double>& v, const std::vector<double>& p, double mq, double dt) {
    std::vector<double> temp = AddVectors(v, DivideVector(p, std::sqrt(mq * mq + DotProduct(p, p))));
    return 30 * dt * std::sqrt(DotProduct(temp, temp)) / 2.0;
}



std::vector<double> CalculateForces(double mq, double Lambda, const std::vector<double>& v, const std::vector<double>& vc, const std::vector<double>& s, double vc0, double vp, int loct, int q, const std::vector<double>& r) {
//total forces on quirks
    
    double beta = 0.0;
    for (double component : v) {
        beta += component * component;
    }
    beta = std::sqrt(beta);
// First term: -Lambda^2 / 100 * sqrt(1 - vc0^2) * s

std::vector<long double> s_long(s.begin(), s.end());
std::vector<long double> vc_long(vc.begin(), vc.end());

    long double vc0_ld = static_cast<long double>(vc0);
long double term1_factor = -Lambda * Lambda / 100.0L * std::sqrt(1.0L - vc0_ld * vc0_ld);
std::vector<long double> term1 = MultiplyVectorLong(s_long, term1_factor);

   // std::vector<double> term1 = MultiplyVector(s, -Lambda * Lambda / 100.0 * std::sqrt(1 - vc0 * vc0));
    //std::cout<<"term1: "<<term1[0]<<", "<<term1[1]<<", "<<term1[2]<<std::endl;
    
long double vp_ld = static_cast<long double>(vp);
long double term2_factor = -Lambda * Lambda / 100.0L * vp_ld / std::sqrt(1.0L - vc0_ld * vc0_ld);
std::vector<long double> term2 = MultiplyVectorLong(vc_long, term2_factor);

// Second term: -Lambda^2 / 100 * vp / sqrt(1 - vc0^2) * vc

//std::vector<double> term2 = MultiplyVector(vc, -Lambda * Lambda / 100.0 * vp / std::sqrt(1 - vc0 * vc0));
//std::cout<<"term2: "<<term2[0]<<", "<<term2[1]<<", "<<term2[2]<<std::endl;
// std::cout<<"vc: "<<vc[0]<<", "<<vc[1]<<", "<<vc[2]<<std::endl;
//std::cout<<"vp: "<<vp<<std::endl;
// Third term: 0.587 * q * Cross(v, Bct(r[0], r[1], r[2]))
// Convert the cross product to std::vector<long double>
std::vector<long double> v_long(v.begin(), v.end());
std::vector<long double> bct_long = BctLong(r[0], r[1], r[2]);

// Use the CrossLong function
std::vector<long double> crossProductLong = CrossLong(v_long, bct_long);

// Multiply the cross product by the scalar
std::vector<long double> term3 = MultiplyVectorLong(crossProductLong, 0.587L * static_cast<long double>(q));


//std::vector<double> term3 = MultiplyVector(Cross(v, Bct(r[0], r[1], r[2])), 0.587 * q);
    //std::cout<<"term3: "<<term3[0]<<", "<<term3[1]<<", "<<term3[2]<<std::endl;
// Fourth term: based on loct value
    double term4;
    switch (loct) {
        case 0:
            term4 =  0.0;
            break;
        case 1:
            term4 = EoxCuAll(mq, 1, beta);
            break;
        case 2:
            term4 = EoxCuAll(mq, 1, beta);
            
            break;
        case 3:
            term4 = EoxCcAll(mq, 1, beta);
            break;
        case 4:
            term4 = EoxRockAll(mq, 1, beta);
            break;
        default:
            throw std::invalid_argument("Invalid loct value");
    }
    //term4 = MultiplyVector(Normalize(v), DotProduct(term4, Normalize(v)));
    std::vector<long double> term4vec = MultiplyVectorLong(NormalizeLong(v), static_cast<long double>(term4));
    //std::cout<<"normalize vector: "<<term4vec[0]<<", "<<term4vec[1]<<", "<<term4vec[2]<<std::endl;
     //   std::cout<<"v : "<<v[0]<<", "<<v[1]<<", "<<v[2]<<std::endl;
   
  
   
    // Summing up all terms
    std::vector<double> force(3);
    for (size_t i = 0; i < force.size(); ++i) {
        force[i] = (static_cast<double>(term1[i]) + static_cast<double>(term2[i]) + static_cast<double>(term3[i]) -  term4vec[i]) / 6.58;
    }

    return force;
}


std::vector<double> CalculateForcesWithGaus(double mq, double Lambda, const std::vector<double>& v, const std::vector<double>& vc, const std::vector<double>& s, double vc0, double vp, int loct, int q, const std::vector<double>& r, double dx) {
//total forces on quirks with gaus de/dex from materials
    double beta = 0.0;
    for (double component : v) {
        beta += component * component;
    }
    beta = std::sqrt(beta);
    // First term: -Lambda^2 / 100 * sqrt(1 - vc0^2) * s
    std::vector<double> term1 = MultiplyVector(s, -Lambda * Lambda / 100.0 * std::sqrt(1 - vc0 * vc0));
    //std::cout<<"term1: "<<term1[0]<<", "<<term1[1]<<", "<<term1[2]<<std::endl;

    // Second term: -Lambda^2 / 100 * vp / sqrt(1 - vc0^2) * vc
    std::vector<double> term2 = MultiplyVector(vc, -Lambda * Lambda / 100.0 * vp / std::sqrt(1 - vc0 * vc0));

//std::cout<<"term2: "<<term2[0]<<", "<<term2[1]<<", "<<term2[2]<<std::endl;
//std::cout<<"vc: "<<vc[0]<<", "<<vc[1]<<", "<<vc[2]<<std::endl;

    // Third term: 0.587 * q * Cross(v, Bct(r[0], r[1], r[2]))
    std::vector<double> term3 = MultiplyVector(Cross(v, Bct(r[0], r[1], r[2])), 0.587 * q);
    //std::cout<<"term3: "<<term3[0]<<", "<<term3[1]<<", "<<term3[2]<<std::endl;
     
    // Fourth term: based on loct value with Gaussian variation
    double term4;
    switch (loct) {
        case 0:
            term4 = 0.0;
            break;
        case 1:
            term4 = EoxCuGaus(mq, 1, beta, dx);
            
            break;
        case 2:
            term4 = EoxCuGaus(mq, 1, beta, dx);
         
            break;
        case 3:
           term4 = EoxCcGaus(mq, 1, beta, dx);
       
            break;
        case 4:
            term4 = EoxRockGaus(mq, 1, beta, dx);
        
            break;
        default:
            throw std::invalid_argument("Invalid loct value");
    }
    //std::cout<<"term4: "<<term4<<std::endl;
    //std::cout<<"beta: "<<beta<<std::endl;
    //term4 = MultiplyVector(Normalize(v), DotProduct(term4, Normalize(v)));
    std::vector<double> term4vec = MultiplyVector(Normalize(v), term4);
 
    // Summing up all terms
    std::vector<double> force(v.size());
    for (size_t i = 0; i < force.size(); ++i) {
        force[i] = (term1[i] + term2[i] + term3[i] - term4vec[i]) / 6.58;
    }

    return force;
}

int main(int argc, char* argv[]) {

    double back = 50e6; // Default back value in micrometers
    double Lambda = 500.0; //Default lambda value in eV
    std::string inputFileName;
    int seed = 0;


    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-b" && i + 1 < argc) {
            back = 1e6*std::atof(argv[++i]); 
        } else if (arg == "-l" && i + 1 < argc) {
            Lambda = std::atof(argv[++i]);
         } else if (arg == "-s" && i + 1 < argc) {
            seed = std::atoi(argv[++i]); 
        } else if (arg[0] != '-') {
            inputFileName = arg; 
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }

     std::mt19937 gen(seed); //initialize RNG with seed

     if (inputFileName.empty()) {
        std::cerr << "Usage: " << argv[0] << " [-b back_value] [-l lambda_value] [-s seed]<input file>" << std::endl;
        return 1;
    }
   
    std::filesystem::path inputPath(inputFileName);
     std::ifstream inputFile(inputFileName);
     std::string stem = inputPath.stem().string();

     // format Lambda without trailing zeros
    std::ostringstream lambdaStream;
    lambdaStream << std::fixed << Lambda;
    std::string lambdaStr = lambdaStream.str();
    lambdaStr.erase(lambdaStr.find_last_not_of('0') + 1);
    if (lambdaStr.back() == '.') {
        lambdaStr.pop_back(); 
    }

    // format back using scientific notation without '+'
    std::ostringstream backStream;
    backStream << std::scientific << back;
    std::string backStr = backStream.str();
    
    //remove the "+" in scientific notation and unnecessary zeros
    size_t ePos = backStr.find('e');
    if (ePos != std::string::npos) {
        std::string exponent = backStr.substr(ePos + 1);
        if (exponent[0] == '+') {
            exponent.erase(0, 1); 
        }
        // construct backStr with significant figures
        std::string significant = backStr.substr(0, ePos);
        significant.erase(significant.find_last_not_of('0') + 1);
        if (significant.back() == '.') {
            significant.pop_back(); 
        }
        backStr = significant + "e" + exponent;
    }

    //output file name for quirks
    std::string outputFileName = stem + "_" + lambdaStr +"eV.txt";
    std::ofstream outputFile(outputFileName);

     //output file name for quirk trajectory
    std::string outputFileTrajectoryName = stem + "_" + lambdaStr +"eV_trajectory.txt";
    std::ofstream outputFileTrajectory(outputFileTrajectoryName );


  

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }


 std::cout << "Running pre-FASER quirk simulation with the following paramters: " << std::endl;
    std::cout << "Final distance: " << back/(1.0e6) <<"m"<< std::endl;
    std::cout << "Lambda: " << Lambda <<"eV"<< std::endl;
   


    // read data from inputFile and perform initial setup
    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(inputFile, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        data.push_back(row);
    }

    for (int h = 0; h <= 0; ++h) {
        auto start = std::chrono::high_resolution_clock::now();
        // Set initial conditions
        double front = 19e6; // in micrometers
        
        int q1 = 1, q2 = -1;
        double mq = round(data[2 * h][6]); // Quirk mass in GeV
        double direc = (data[2 * h][4] + data[2 * h + 1][4] > 0) ? 1 : -1;

        

        std::vector<double> p1 = {direc * data[2 * h ][2], direc * data[2 * h ][3], direc * data[2 * h ][4]};
        std::vector<double> p2 = {direc * data[2 * h + 1][2], direc * data[2 * h + 1][3], direc * data[2 * h + 1][4]};

        
        double E1 = sqrt(mq * mq + p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
        double E2 = sqrt(mq * mq + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);

        std::vector<double> v1 = {p1[0] / E1, p1[1] / E1, p1[2] / E1};
        std::vector<double> v2 = {p2[0] / E2, p2[1] / E2, p2[2] / E2};

        double v10 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        double v20 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);

        std::vector<double> Beta = {(p1[0] + p2[0]) / (E1 + E2), (p1[1] + p2[1]) / (E1 + E2), (p1[2] + p2[2]) / (E1 + E2)};

        double t1q = 658 * ((2 * mq) / (Lambda * Lambda)) * sqrt(pow((E1 + E2) / (2 * mq), 2) - 1 / (1 - (Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2])));
        double dt = std::min(0.03, t1q / 10000);

        int nsf = floor(front / (3e5 * t1q * Beta[2]));
        int ns = (nsf % 2 == 0) ? nsf : nsf - 1;
        double t1 = ns * t1q;
        double t2 = t1;

        std::vector<double> r1 = {3e5 * ns * t1q * Beta[0], 3e5 * ns * t1q * Beta[1], 3e5 * ns * t1q * Beta[2]};
        std::vector<double> r2 = r1;

        int stepcount = 0;
        int n = 1;
        // while (!(((r1[2] > back) && (r2[2] > back)) 

        double lastSaveTime = 0;
        double saveInterval = .1; 

        double dx1pre = 0.0, dx2pre = 0.0; 

        while (!( (sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2]) < 0.01))) {
            int loct1 = Loct(r1[0], r1[1], r1[2]);
            int loct2 = Loct(r2[0], r2[1], r2[2]);

           

            int layer1i = Layer(r1[0], r1[1], r1[2]);
            int layer2i = Layer(r2[0], r2[1], r2[2]);

            stepcount++;
            

            // Recalculate energies and velocities
            E1 = sqrt(mq * mq + p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
            v1 = {p1[0] / E1, p1[1] / E1, p1[2] / E1};
            v10 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

            

            E2 = sqrt(mq * mq + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);
            v2 = {p2[0] / E2, p2[1] / E2, p2[2] / E2};
            v20 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);



            Beta = {(p1[0] + p2[0]) / (E1 + E2), (p1[1] + p2[1]) / (E1 + E2), (p1[2] + p2[2]) / (E1 + E2)};

            // Direction of the string at quirk 1
            std::vector<double> s1(3), s2(3);
            if (n == 1) {
                s1 = Normalize(SubtractVectors(v1, v2));
            } else {
                s1 = Normalize(AddVectors(MultiplyVector(SubtractVectors(SubtractVectors(r1, r2), MultiplyVector(Beta, 300000 * (t1 - t2))), (1 - DotProduct(Beta, Beta))), MultiplyVector(SubtractVectors(Beta, v1), DotProduct(SubtractVectors(SubtractVectors(r1, r2), MultiplyVector(Beta, 300000 * (t1 - t2))), Beta))));
            }

            double vp1 = DotProduct(v1, s1);
            std::vector<double> vc1 = SubtractVectors(v1, MultiplyVector(s1, vp1));
            double vc10 = sqrt(DotProduct(vc1, vc1));

            // Direction of the string at quirk 2
            if (n == 1) {
                s2 = Normalize(SubtractVectors(v2, v1));
            } else {
                s2 =  Normalize(AddVectors(MultiplyVector(SubtractVectors(SubtractVectors(r2, r1), MultiplyVector(Beta, 300000 * (t2 - t1))), (1 - DotProduct(Beta, Beta))), MultiplyVector(SubtractVectors(Beta, v2), DotProduct(SubtractVectors(SubtractVectors(r2, r1), MultiplyVector(Beta, 300000 * (t2 - t1))), Beta))));
            }

            double vp2 = DotProduct(v2, s2);
            std::vector<double> vc2 = SubtractVectors(v2, MultiplyVector(s2, vp2));
            double vc20 = sqrt(DotProduct(vc2, vc2));

            // Estimation of the travel distance (in cm) using the average de/dx value
           
            if (loct1 > 0 || loct2 > 0) {

            
                std::vector<double> F1pre = CalculateForces(mq, Lambda, v1, vc1, s1, vc10, vp1, loct1, q1, r1);
                std::vector<double> F2pre = CalculateForces(mq, Lambda, v2, vc2, s2, vc20, vp2, loct2, q2, r2);

             

                double ct1pre = CalculateCt(v1, Beta, r1, r2, F1pre, E1, E2);
                double ct2pre = CalculateCt(v2, Beta, r2, r1, F2pre, E1, E2);

                double dt1pre, dt2pre;
                if (abs(ct1pre) < abs(ct2pre)) {
                    dt1pre = dt;
                    dt2pre = dt1pre * ct1pre / ct2pre;
                } else {
                    dt2pre = dt;
                    dt1pre = dt2pre * ct2pre / ct1pre;
                }


                std::vector<double> p1pre = AddVectors(p1, MultiplyVector(F1pre, dt1pre));
                std::vector<double> p2pre = AddVectors(p2, MultiplyVector(F2pre, dt2pre));

                dx1pre = CalculateDistance(v1, p1pre, mq, dt1pre);
                dx2pre = CalculateDistance(v2, p2pre, mq, dt2pre);
            }
         
            // Recalculate forces using normally distributed de/dx
            std::vector<double> F1 = CalculateForcesWithGaus(mq, Lambda, v1, vc1, s1, vc10, vp1, loct1, q1, r1, dx1pre);
            std::vector<double> F2 = CalculateForcesWithGaus(mq, Lambda, v2, vc2, s2, vc20, vp2, loct2, q2, r2, dx2pre);
       

            double ct1 = CalculateCt(v1, Beta, r1, r2, F1, E1, E2);
            double ct2 = CalculateCt(v2, Beta, r2, r1, F2, E1, E2);

            double dt1, dt2;
     
            if (abs(ct1) < abs(ct2)) {
                dt1 = dt;
                dt2 = dt1 * ct1 / ct2;
            } else {
                dt2 = dt;
                dt1 = dt2 * ct2 / ct1;
            
            }

            // Update quirk momentum and position
            p1 = AddVectors(p1, MultiplyVector(F1, dt1));
            p2 = AddVectors(p2, MultiplyVector(F2, dt2));

            r1 = AddVectors(r1, MultiplyVector(AddVectors(v1, DivideVector(p1, sqrt(mq * mq + DotProduct(p1, p1)))), 300000.0 * dt1 / 2));
            r2 = AddVectors(r2, MultiplyVector(AddVectors(v2, DivideVector(p2, sqrt(mq * mq + DotProduct(p2, p2)))), 300000.0 * dt2 / 2));

            t1 += dt1;
            t2 += dt2;

            // Determine the detector scintillators
            int layer1f = Layer(r1[0], r1[1], r1[2]);
            int layer2f = Layer(r2[0], r2[1], r2[2]);

        
        //currently breaking here to locate start and source of numerical instability
         /*if (n==253056){
            std::cout<<std::setprecision(16);
            std::cout<<"F1: "<<F1[0]<<", "<<F1[1]<<", "<<F1[2]<<std::endl;
            std::cout<<"F2: "<<F2[0]<<", "<<F2[1]<<", "<<F2[2]<<std::endl;
            
            std::cout<<"r1: "<<r1[0]<<", "<<r1[1]<<", "<<r1[2]<<std::endl;
            std::cout<<"r2: "<<r2[0]<<", "<<r2[1]<<", "<<r2[2]<<std::endl;
            std::cout<<"E1: "<<E1<<std::endl;
            std::cout<<"E1: "<<E2<<std::endl;
            std::cout<<"p1: "<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<std::endl;
            std::cout<<"p2: "<<p2[0]<<", "<<p2[1]<<", "<<p2[2]<<std::endl;

            std::cout<<"s1: "<<s1[0]<<", "<<s1[1]<<", "<<s1[2]<<std::endl;
            std::cout<<"s2: "<<s2[0]<<", "<<s2[1]<<", "<<s2[2]<<std::endl;
    
            std::cout<<"n: "<<n<<std::endl;

            std::cout<<"dt1: "<<dt1<<std::endl;
            std::cout<<"dt2: "<<dt2<<std::endl;

            std::cout<<"t1: "<<t1<<std::endl;
            std::cout<<"t2: "<<t2<<std::endl;

            std::cout<<"ct1: "<<ct1<<std::endl;
            std::cout<<"ct2: "<<ct2<<std::endl;


            
            std::cout<<"loct1: "<<loct1<<std::endl;
            std::cout<<"loct2: "<<loct2<<std::endl;
            std::cout<<"vp2: "<<vp2<<std::endl;
            std::cout<<"dx2pre: "<<dx2pre<<std::endl;
            std::cout<<"vc20: "<<vc20<<std::endl;
            std::cout<<"s2: "<<s2[0]<<", "<<s2[1]<<", "<<s2[2]<<std::endl;
            std::cout<<"vc2: "<< vc2[0]<<", "<<vc2[1]<<", "<<vc2[2]<<std::endl;
            std::cout<<"v2: "<< v2[0]<<", "<<v2[1]<<", "<<v2[2]<<std::endl;
            std::cout<<"beta: "<<sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2])<<std::endl;

             break;
            }*/
        //save trajectory info
        if(t1 - lastSaveTime >= saveInterval){
            outputFileTrajectory << std::setprecision(16) <<  t1 << " " << r1[0] << " " << r1[1] << " " << r1[2] << " " << "\n";
            outputFileTrajectory << std::setprecision(16) <<  t1 << " " << r2[0] << " " << r2[1] << " " << r2[2] << " " << "\n";
            lastSaveTime = t1;
        }

     if (stepcount%1000000==0){
  
        std::cout<<"z1: "<<r1[2]<<std::endl;
         std::cout<<"beta: "<<sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2])<<std::endl;

     }


    if (r1[2] > back) {

            // output hit information if condition met
  
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;

        outputFile << std::setprecision(16) << h << " " << mq << " " << Lambda << " " << t1q << " 1 " << 1 << " " << t1 << " "
                   << r1[0] << " " << r1[1] << " " << r1[2] << " "<< p1[0] << " " << p1[1] << " " << p1[2] <<" " <<duration.count()  << "\n";
        outputFile << std::setprecision(16) << h << " " << mq << " " << Lambda << " " << t1q << " 2 " << 1 << " " << t2 << " "
                   << r2[0] << " " << r2[1] << " " << r2[2] << " "<< p2[0] << " " << p2[1] << " " << p2[2] <<" " << duration.count() << "\n";

      
       break;
    }
     n++;

   
}
    std::cout<<h<<std::endl;
}
 

    inputFile.close();
    outputFile.close();
    outputFileTrajectory.close();

    return 0;
}