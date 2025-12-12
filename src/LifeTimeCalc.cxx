
#include "../include/LifeTimeCalc.h"

void FourVector::boost(const std::array<double,3>& beta) {
    double b2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bp = beta[0]*px + beta[1]*py + beta[2]*pz;
    double factor = (gamma - 1.0) * bp / b2 - gamma * e;
    px += factor * beta[0];
    py += factor * beta[1];
    pz += factor * beta[2];
    e  = gamma * (e - bp);
}

FourVector Boost_lab_to_com(double p1x,double p1y,double p1z,double e1,
                            double p2x,double p2y,double p2z,double e2)
{
    FourVector v1(p1x,p1y,p1z,e1);
    FourVector v2(p2x,p2y,p2z,e2);
    FourVector tot = v1 + v2;
    auto b = tot.boost_vector();
    v1.boost({-b[0], -b[1], -b[2]});
    return v1;
}

//decay distance
double DecayDistance(double Epsilon, double Epsilon1, double Lambda_eV, double m, const std::string& name,
                    double p1x, double p1y, double p1z,
                    double p2x, double p2y, double p2z)
{
    const double hbar = 6.5821e-25; // GeV*s
    const double c = 2.99792458e8;  // m/s
    const double alpha_EM = 1.0 / 137.0;
    const double Nic = 2.0;
    const double eV_to_GeV = 1e-9;
    double Lambda_IC = Lambda_eV * eV_to_GeV; 

    double e1 = std::sqrt(p1x*p1x + p1y*p1y + p1z*p1z + m*m);
    double e2 = std::sqrt(p2x*p2x + p2y*p2y + p2z*p2z + m*m);
    FourVector com = Boost_lab_to_com(p1x,p1y,p1z,e1,p2x,p2y,p2z,e2);

    double plx = p1x + p2x;
    double ply = p1y + p2y;
    double plz = p1z + p2z;
    double kl = std::sqrt(plx*plx + ply*ply + plz*plz);
    double kc = std::sqrt(com.px*com.px + com.py*com.py + com.pz*com.pz);

    double el = std::sqrt(kl*kl + (2*m)*(2*m));
    double ec = std::sqrt(kc*kc + m*m);
    double v = std::fabs(kl / el) * c;

    // quirk models
    int ns = 0, nf = 1;
    double Qt = 1.0, NC = 1.0;
    if (name == "f31") { Qt = 1.0; }
    else if (name == "s31") { ns=1; nf=0; Qt=1.0; }
    else if (name == "f33") { ns=0; nf=1; NC=3; Qt=2.0/3.0; }
    else if (name == "s33") { ns=1; nf=0; NC=3; Qt=2.0/3.0; }

    double b0 = (1.0/(4*M_PI))*((11.0/3.0)*Nic - (2.0/3.0)*nf - (1.0/6.0)*ns);
    double alpha_IC = 1.0 / (b0 * std::log(m*m / (Lambda_IC*Lambda_IC)));

    double EL = 2*ec - 2*m;
    double T_gl = (4*std::sqrt(2*m)*((std::pow(EL,1.5)-std::pow(std::sqrt(m*Lambda_IC),1.5))
                 /(3*std::pow(Lambda_IC,3)*Epsilon))) * hbar;
    double T_gc = (2*M_PI*alpha_IC*(std::sqrt(m)/(std::pow(Lambda_IC,1.5)*Epsilon))
                 - (2*M_PI/(Epsilon*Lambda_IC))) * hbar;

    double T_el = ((3*(EL - std::sqrt(m*Lambda_IC))*m*m) /
                  (16*Qt*Qt*M_PI*alpha_EM*std::pow(Lambda_IC,4))) * hbar;
    double T_ec = ((alpha_IC*alpha_IC*m*m) /
                  (16*Qt*Qt*M_PI*alpha_EM*std::pow(Lambda_IC,3))
                  - (1.0 / (16*Qt*Qt*M_PI*alpha_EM*std::pow(alpha_IC,4)*m))) * hbar;

    double t00 = (1.0 / (1.0/T_gl + 1.0/T_el)) + (1.0 / (1.0/T_gc + 1.0/T_ec));
    double tau_lab = t00 / std::sqrt(1 - (v*v)/(c*c));

    // decay distance
    double decay_distance = v * tau_lab;
    return decay_distance;
}

//standard deviation based on lab frame 4-momenta and lambda
double DecayStandardDeviation(double Lambda_eV, double m,
                         double p1x, double p1y, double p1z, double e1,
                         double p2x, double p2y, double p2z, double e2)
{
    const double c = 2.99792458e8;  // m/s
    const double eV_to_GeV = 1e-9;
    double Lambda_IC = Lambda_eV * eV_to_GeV;  // Convert to GeV
    
    // Boost to CoM frame
    FourVector com = Boost_lab_to_com(p1x, p1y, p1z, e1, p2x, p2y, p2z, e2);
    
    // Calculate CoM frame quantities
    double kc = std::sqrt(com.px*com.px + com.py*com.py + com.pz*com.pz);
    double ec = std::sqrt(kc*kc + m*m);
    double EL = 2*ec - 2*m;  // Total kinetic energy in CoM frame
    
    // Calculate lab frame quantities for beta
    double plx = p1x + p2x;
    double ply = p1y + p2y;
    double plz = p1z + p2z;
    double kl = std::sqrt(plx*plx + ply*ply + plz*plz);
    double el = std::sqrt(kl*kl + (2*m)*(2*m));
    double v = std::fabs(kl / el) * c;
    double beta_cm = v / c;
    double beta_sq = beta_cm * beta_cm;
    
    
    // Calculate oscillation period (t1q)
    // t1q = 658 * (2*m / Lambda^2) * sqrt(((E1+E2)/(2*m))^2 - 1/(1-beta^2))
    // Lambda is in eV, result is in nanoseconds
    double E_total = e1 + e2;
    double E_over_2m = E_total / (2 * m);
    double gamma_term = E_over_2m*E_over_2m - 1.0 / (1.0 - beta_sq);
    if (gamma_term < 0) {
        gamma_term = 0;  
    }
    double gamma_factor = std::sqrt(gamma_term);
    
    // t1q in nanoseconds, convert to seconds
    // Lambda_eV is already in eV
    double t1q_ns = 658.0 * ((2.0 * m) / (Lambda_eV * Lambda_eV)) * gamma_factor;
    double t1q_seconds = t1q_ns * 1e-9;  // Convert nanoseconds to seconds
    
    // Calculate standard deviation: sqrt(N) * oscillation_period
    // N = (total kinetic energy in CoM frame) / Lambda
    double N = EL / Lambda_IC;
    double std_dev_time = std::sqrt(N) * t1q_seconds;
    
    // Convert to distance 
    double std_dev_distance = v * std_dev_time;
    
    return std_dev_distance;
}



// All args in meters
// returns P(L > Lcut) for L ~ N(mean, sigma^2)
 double LifetimeSurvivalProbGauss(double mean,double sigma,double Lcut)
{
    // Bad sigma -> treat as step at mean
    if (sigma <= 0.0 || !std::isfinite(sigma)) {
        return (mean > Lcut) ? 1.0 : 0.0;
    }

    const double z   = (Lcut - mean) / sigma;
    const double arg = z / std::sqrt(2.0);

    // P(L > Lcut) = 0.5 * erfc( (Lcut - mean) / (sqrt(2)*sigma) )
    double p = 0.5 * std::erfc(arg);

    if (!std::isfinite(p)) p = 0.0;
    p = std::max(0.0, std::min(1.0, p));
    return p;
}
