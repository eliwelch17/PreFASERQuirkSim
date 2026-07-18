#include "../include/QuirkDynamics.h"

#include "../include/DetectorGeometry.h"
#include "../include/EnergyLoss.h"
#include "../include/MagneticField.h"
#include "../include/VectorUtils.h"

#include <cmath>
#include <stdexcept>

double CalculateCt(const std::vector<double> &v, const std::vector<double> &Beta, const std::vector<double> &r1, const std::vector<double> &r2, const std::vector<double> &F, double E1, double E2)
{
    return 1 - DotProduct(v, Beta) - DotProduct(SubtractVectors(r1, r2), SubtractVectors(F, MultiplyVector(Beta, DotProduct(v, F)))) / (300000 * (E1 + E2));
}

// Function to calculate the travel distance
double CalculateDistance(const std::vector<double> &v, const std::vector<double> &p, double mq, double dt)
{
    std::vector<double> temp = AddVectors(v, DivideVector(p, std::sqrt(mq * mq + DotProduct(p, p))));
    return 30 * dt * std::sqrt(DotProduct(temp, temp)) / 2.0;
}

std::vector<double> CalculateForces(double mq, double Lambda, const std::vector<double> &v, const std::vector<double> &vc, const std::vector<double> &s, double vc0, double vp, int loct, int q, const std::vector<double> &r)
{
    // total forces on quirks

    double beta = 0.0;
    for (double component : v)
    {
        beta += component * component;
    }
    beta = std::sqrt(beta);
    // First term: -Lambda^2 / 100 * sqrt(1 - vc0^2) * s


    double term1_factor = -Lambda * Lambda / 100.0L * std::sqrt(1.0L - vc0* vc0);
    std::vector<double> term1 = MultiplyVector(s, term1_factor);

    // std::vector<double> term1 = MultiplyVector(s, -Lambda * Lambda / 100.0 * std::sqrt(1 - vc0 * vc0));
    // std::cout<<"term1: "<<term1[0]<<", "<<term1[1]<<", "<<term1[2]<<std::endl;

    double term2_factor = -Lambda * Lambda / 100.0L * vp / std::sqrt(1.0L - vc0 * vc0);
    std::vector<double> term2 = MultiplyVector(vc, term2_factor);

    // Second term: -Lambda^2 / 100 * vp / sqrt(1 - vc0^2) * vc

 
    std::vector<double> bct = BctAdv(all_magnets, r[0], r[1], r[2]);
    std::vector< double> crossProduct = Cross(v, bct);
    std::vector<double> term3 = MultiplyVector(crossProduct, 0.587L * q);
     
    // Fourth term: based on loct value
    double term4;
    switch (loct)
    {
    case 0:
        term4 = 0.0;
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

    std::vector<double> term4vec = MultiplyVector(Normalize(v), term4);

    // Summing up all terms
    std::vector<double> force(3);
    for (size_t i = 0; i < force.size(); ++i)
    {
        force[i] = (static_cast<double>(term1[i]) + static_cast<double>(term2[i]) + static_cast<double>(term3[i]) - term4vec[i]) / 6.58;
    }

    return force;
}

std::vector<double> CalculateForcesWithGaus(double mq, double Lambda, const std::vector<double> &v, const std::vector<double> &vc, const std::vector<double> &s, double vc0, double vp, int loct, int q, const std::vector<double> &r, double dx, std::mt19937 &gen)
{
    // total forces on quirks with gaus de/dex from materials
    
    double beta = 0.0;
    for (double component : v)
    {
        beta += component * component;
    }
    beta = std::sqrt(beta);
    // First term: -Lambda^2 / 100 * sqrt(1 - vc0^2) * s
    double sqrtTerm1 = std::sqrt(static_cast<double>(1.0) - vc0 * vc0);
    std::vector<double> term1 = MultiplyVector(s, -Lambda * Lambda / 100.0 * static_cast<double>(sqrtTerm1));

    // Second term: -Lambda^2 / 100 * vp / sqrt(1 - vc0^2) * vc
    double sqrtTerm2 = std::sqrt(1.0 - vc0 * vc0);
    std::vector<double> term2 = MultiplyVector(vc, -Lambda * Lambda / 100.0 * vp / sqrtTerm2);

    // Third term: 0.587 * q * Cross(v, Bct(r[0], r[1], r[2]))
    
    std::vector<double> bct = BctAdv(all_magnets, r[0], r[1], r[2]);
    // std::vector<long double> bct_long = BctLong(r[0], r[1], r[2]);

   
    std::vector<double> term3 = MultiplyVector(Cross(v, bct), 0.587L * q);

    // Fourth term: based on loct value with Gaussian variation
    double term4;
    switch (loct)
    {
    case 0:
        term4 = 0.0;
        break;
    case 1:
        term4 = EoxCuGaus(mq, 1, beta, dx, gen);

        break;
    case 2:
        term4 = EoxCuGaus(mq, 1, beta, dx, gen);

        break;
    case 3:
        term4 = EoxCcGaus(mq, 1, beta, dx, gen);

        break;
    case 4:
        term4 = EoxRockGaus(mq, 1, beta, dx, gen);

        break;
    default:
        throw std::invalid_argument("Invalid loct value");
    }

    std::vector<double> term4vec = MultiplyVector(Normalize(v), term4);
    // Summing up all terms
    std::vector<double> force(3);
    for (size_t i = 0; i < force.size(); ++i)
    {
        force[i] = (term1[i] + term2[i] + term3[i] - term4vec[i]) / 6.58;
    }

    return force;
}




void SyncQuirksToSameTime(
    double t_star,
    std::vector<double>& r1, std::vector<double>& p1, double& t1, int q1,
    std::vector<double>& r2, std::vector<double>& p2, double& t2, int q2,
    double mq, double Lambda)
{
  auto v_from_p = [&](const std::vector<double>& p){
    double E = std::sqrt(mq*mq + DotProduct(p,p));
    return DivideVector(p, E);
  };

  auto advance_lag = [&](std::vector<double>& r,
                         std::vector<double>& p,
                         double& t,
                         const std::vector<double>& rO,
                         const std::vector<double>& pO,
                         int q,
                         double dt)
  {
    if (dt <= 0) return;

    // velocity at start
    std::vector<double> v = v_from_p(p);

    // string direction: unit vector from other quirk to this quirk
    std::vector<double> dr = SubtractVectors(r, rO);
    std::vector<double> s  = Normalize(dr);

    // decompose v into parallel / perpendicular to s
    double vpar = DotProduct(v, s);
    std::vector<double> vperp = SubtractVectors(v, MultiplyVector(s, vpar));
    double vperp0 = std::sqrt(DotProduct(vperp, vperp));

    // location/material
    int loct = Loct(r[0], r[1], r[2]);

   
    std::vector<double> F = CalculateForces(mq, Lambda, v, vperp, s, vperp0, vpar, loct, q, r);

    // trapezoidal update over dt
    std::vector<double> p_new = AddVectors(p, MultiplyVector(F, dt));
    std::vector<double> v_new = v_from_p(p_new);
    r = AddVectors(r, MultiplyVector(AddVectors(v, v_new), 300000.0 * dt * 0.5));
    p = p_new;
    t += dt;
  };

  if (t1 < t_star) advance_lag(r1, p1, t1, r2, p2, q1, t_star - t1);
  if (t2 < t_star) advance_lag(r2, p2, t2, r1, p1, q2, t_star - t2);
}
