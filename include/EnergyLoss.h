#pragma once

#include <functional>
#include <random>
#include <vector>

double Gamma(double beta);
double Emax(double m, double beta);
double Xi(double z, double ZoA, double rho, double beta);
double Delta(double x, double x0, double x1, double Cbar, double a, double k, double d0);
double Eox(double m, double z, double ZoA, double rho, double I0, double x0, double x1,
           double Cbar, double a, double k, double d0, double beta);
double EoxCuAll(double mq, int param, double beta);
double EoxCcAll(double mq, int param, double beta);
double EoxRockAll(double mq, int param, double beta);
double EoxGaus(double m, double z, double ZoA, double rho, double I0, double x0, double x1,
               double Cbar, double a, double k, double d0, double beta, double delta_x,
               std::function<double(double, int, double)> EoxAllFunc, std::mt19937& gen);
double EoxCuGaus(double mq, int param, double beta, double dx, std::mt19937& gen);
double EoxCcGaus(double mq, int param, double beta, double dx, std::mt19937& gen);
double EoxRockGaus(double mq, int param, double beta, double dx, std::mt19937& gen);
std::vector<double> EoxCu(double mq, int param, const std::vector<double>& v);
std::vector<double> EoxCc(double mq, int param, const std::vector<double>& v);
std::vector<double> EoxRock(double mq, int param, const std::vector<double>& v);
