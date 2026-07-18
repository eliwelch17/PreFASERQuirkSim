#include "../include/VectorUtils.h"

#include <cmath>
#include <stdexcept>

std::vector<double> Cross(const std::vector<double> &v1, const std::vector<double> &v2)
{
    // cross product function
    if (v1.size() != 3 || v2.size() != 3)
    {
        throw std::invalid_argument("Both vectors must be 3-dimensional");
    }
    std::vector<double> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

std::vector<long double> CrossLong(const std::vector<long double> &v1, const std::vector<long double> &v2)
{
    // cross product function with long's
    if (v1.size() != 3 || v2.size() != 3)
    {
        throw std::invalid_argument("Both vectors must be 3-dimensional");
    }
    std::vector<long double> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

std::vector<double> SubtractVectors(const std::vector<double> &v1, const std::vector<double> &v2)
{
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

std::vector<double> AddVectors(const std::vector<double> &v1, const std::vector<double> &v2)
{
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

std::vector<double> MultiplyVector(const std::vector<double> &v, double scalar)
{
    return {v[0] * scalar, v[1] * scalar, v[2] * scalar};
}



std::vector<double> DivideVector(const std::vector<double>& v, double s){
  if (!std::isfinite(s) || std::abs(s) < 1e-15) s = (s>=0?1e-15:-1e-15);
  return {v[0]/s, v[1]/s, v[2]/s};
}


double DotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
    double sum=0.0, c=0.0;
    for (size_t i=0;i<v1.size();++i) {
        double y = v1[i]*v2[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

std::vector<double> Normalize(const std::vector<double>& v){
  double n = std::sqrt(DotProduct(v,v));
  if (n < 1e-15) n = 1e-15;
  return {v[0]/n, v[1]/n, v[2]/n};
}


std::vector<double> BoostToCOM(const std::vector<double>& p, double E, const std::vector<double>& Beta) {
    double beta_sq = DotProduct(Beta, Beta);
    if (beta_sq < 1e-15) {
        return p;
    }
    double gamma = 1.0 / std::sqrt(1.0 - beta_sq);
    double p_dot_beta = DotProduct(p, Beta);
    std::vector<double> p_COM = AddVectors(p, MultiplyVector(Beta, (gamma - 1.0) * p_dot_beta / beta_sq - gamma * E));
    return p_COM;
}
