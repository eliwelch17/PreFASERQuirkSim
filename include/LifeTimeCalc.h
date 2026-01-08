#ifndef LIFETIME_CALC_H
#define LIFETIME_CALC_H

#include <cmath>
#include <string>
#include <array>

struct FourVector {
    double px, py, pz, e;
    FourVector(double x=0, double y=0, double z=0, double E=0) : px(x), py(y), pz(z), e(E) {}

    FourVector operator+(const FourVector& o) const {
        return FourVector(px + o.px, py + o.py, pz + o.pz, e + o.e);
    }

    std::array<double,3> boost_vector() const {
        return {-px / e, -py / e, -pz / e};
    }

    void boost(const std::array<double,3>& beta);
};

FourVector Boost_lab_to_com(double p1x, double p1y, double p1z, double e1,double p2x, double p2y, double p2z, double e2);

double DecayDistance(double Epsilon, double Epsilon1, double Lambda_eV, double m, const std::string& name,double p1x, double p1y, double p1z,double p2x, double p2y, double p2z);

double DecayStandardDeviation(double Lambda_eV, double m,double p1x, double p1y, double p1z, double e1,double p2x, double p2y, double p2z, double e2);

// RMS transverse deflection width at tracker plane for quirk-pair system, in micrometers
// sigma_R(epsilon) = d * 0.34 * sqrt(epsilon * N_osc) * (Lambda / |p_tot|)
// where N_osc = EL / Lambda (EL is total CoM kinetic energy), Lambda in eV input is converted internally.
double SigmaR_um(double epsilon, double d_m, double Lambda_eV, double m,
                 double p1x, double p1y, double p1z,
                 double p2x, double p2y, double p2z);


double LifetimeSurvivalProbGauss(double mean,double sigma,double Lcut);
#endif // LIFETIME_CALC_H
