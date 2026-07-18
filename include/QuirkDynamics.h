#pragma once

#include <random>
#include <vector>

double CalculateCt(const std::vector<double>& v, const std::vector<double>& Beta,
                   const std::vector<double>& r1, const std::vector<double>& r2,
                   const std::vector<double>& F, double E1, double E2);
double CalculateDistance(const std::vector<double>& v, const std::vector<double>& p,
                         double mq, double dt);
std::vector<double> CalculateForces(double mq, double Lambda, const std::vector<double>& v,
                                    const std::vector<double>& vc, const std::vector<double>& s,
                                    double vc0, double vp, int loct, int q,
                                    const std::vector<double>& r);
std::vector<double> CalculateForcesWithGaus(double mq, double Lambda,
                                            const std::vector<double>& v,
                                            const std::vector<double>& vc,
                                            const std::vector<double>& s,
                                            double vc0, double vp, int loct, int q,
                                            const std::vector<double>& r, double dx,
                                            std::mt19937& gen);
void SyncQuirksToSameTime(double t_star,
                          std::vector<double>& r1, std::vector<double>& p1, double& t1, int q1,
                          std::vector<double>& r2, std::vector<double>& p2, double& t2, int q2,
                          double mq, double Lambda);
