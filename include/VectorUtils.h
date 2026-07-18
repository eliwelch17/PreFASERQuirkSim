#pragma once

#include <vector>

std::vector<double> Cross(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<long double> CrossLong(const std::vector<long double>& v1,
                                   const std::vector<long double>& v2);
std::vector<double> SubtractVectors(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> AddVectors(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> MultiplyVector(const std::vector<double>& v, double scalar);
std::vector<double> DivideVector(const std::vector<double>& v, double s);
double DotProduct(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> Normalize(const std::vector<double>& v);
std::vector<double> BoostToCOM(const std::vector<double>& p, double E,
                               const std::vector<double>& Beta);
