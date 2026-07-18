#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <deque>
#include <functional>
#include <iomanip>
#include <string>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <limits>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <array>
#include <map>

#include "include/DetectorGeometry.h"
#include "include/LifeTimeCalc.h"
#include "include/MagneticField.h"
#include "include/QuirkDynamics.h"
#include "include/RangeTableLoss.h"
#include "include/VectorUtils.h"

using namespace std;

// ##############################################################################################################################
//
//
//
// PreFaser Quirk Simulation
//
//
//
// ##############################################################################################################################

int main(int argc, char *argv[])
{

    initializeFieldMaps();
   
    double back = 474.6e6-0.2e6; //Vetonu scintillatorin faser coords-.2m buffer. Default back value in micrometers (vetoNu0), faser z=0 at 477.76
    double Lambda = 500.0; // Default lambda value in eV
    std::string inputFileName;
    int seed = 0;        // random seed
    int nquirks = -1;    // number of quirks to simulate
    int divStep = 10000; // not number of steps, step size is proportioal to Lambda squared/stepDen
    bool traj = false;   // trajectory output flag
    int skip = 0;
    int runNum = 0;
    double front = 19e6; // in micrometers
    double beta_cutOff = 0.1; // minimum beta to continue simulating
 

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-b" && i + 1 < argc)
        {
            back = 1e6 * std::atof(argv[++i]);
        }
        else if (arg == "-l" && i + 1 < argc)
        {
            Lambda = std::atof(argv[++i]);
        }
        else if (arg == "-betaCut" && i + 1 < argc)
        {
            beta_cutOff = std::atof(argv[++i]);
        }
        else if (arg == "-skip" && i + 1 < argc)
        {
            skip = std::atof(argv[++i]);
        }
        else if (arg == "-r" && i + 1 < argc)
        {
            runNum = std::atof(argv[++i]);
        }
        else if (arg == "-n" && i + 1 < argc)
        {
            nquirks = std::atof(argv[++i]);
        }
        else if (arg == "-f" && i + 1 < argc)
        {
            front = 1e6 * std::atof(argv[++i]); //input in meters
        }
        else if (arg == "-s" && i + 1 < argc)
        {
            seed = std::atoi(argv[++i]);
        }
        else if (arg == "-t" && i + 1 < argc)
        {
            traj = true;
        }
        else if (arg == "-d" && i + 1 < argc)
        {
            divStep = std::atoi(argv[++i]);
        }
        else if (arg[0] != '-')
        {
            inputFileName = arg;
        }
        else
        {
            std::cerr << "Unknown option: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0] << " [-b <back_value>] [-betaCut <min Beta cut off>] [-l <lambda_value>] [-s <seed>] [-f <front>] [-n <# quirks>] [-d <stepsize divider>] [-skip <# events to skip>] [-runNum <run number>] [-t (trajectory output flag)] <input file>" << std::endl;
            return 1;
        }
    }

    std::mt19937 gen(seed); // initialize RNG with seed

    if (inputFileName.empty())
    {
        std::cerr << "Usage: " << argv[0] << " [-b <back_value>] [-betaCut <min Beta cut off>] [-l <lambda_value>] [-s <seed>] [-n <# quirks>] [-d <stepsize divider>] [-f <front>] [-skip <# events to skip>] [-runNum <run number>] [-t (trajectory output flag)] <input file>" << std::endl;
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
    if (lambdaStr.back() == '.')
    {
        lambdaStr.pop_back();
    }

 
    std::ostringstream backStream;
    backStream << std::scientific << back;
    std::string backStr = backStream.str();


    size_t ePos = backStr.find('e');
    if (ePos != std::string::npos)
    {
        std::string exponent = backStr.substr(ePos + 1);
        if (exponent[0] == '+')
        {
            exponent.erase(0, 1);
        }

        std::string significant = backStr.substr(0, ePos);
        significant.erase(significant.find_last_not_of('0') + 1);
        if (significant.back() == '.')
        {
            significant.pop_back();
        }
        backStr = significant + "e" + exponent;
    }

    // output file name for quirks
    std::string outputFileName = stem + "_" + lambdaStr + "eV" + "_sd_" + std::to_string(seed) + "_" + std::to_string(runNum) + ".txt";
    std::ofstream outputFile(outputFileName);

    std::unique_ptr<std::ofstream> outputFileTrajectory;
    if (traj)
    {
        // output file name for quirk trajectory
        std::string outputFileTrajectoryName = stem + "_" + lambdaStr + "eV_trajectory.txt";
        outputFileTrajectory = std::make_unique<std::ofstream>(outputFileTrajectoryName);
    }

    if (!inputFile.is_open() || !outputFile.is_open())
    {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }

    std::cout << "*****************************************************" << std::endl;
    std::cout << "Running pre-FASER quirk simulation with the following paramters: " << std::endl;
    std::cout << "Final distance: " << back / (1.0e6) << "m" << std::endl;
    std::cout << "Lambda: " << Lambda << "eV" << std::endl;
    std::cout << "Beta cut-off: " << beta_cutOff  << std::endl;
    if (nquirks == -1)
    {
        std::cout << "Number of quirks: All" << std::endl;
    }
    else
    {
        std::cout << "Number of quirks: " << nquirks << std::endl;
    }
    std::cout << "*****************************************************" << std::endl;

    // read data from inputFile and perform initial setup
    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(inputFile, line))
    {
        std::vector<double> row;
        std::stringstream ss(line);
        double value;
        while (ss >> value)
        {
            row.push_back(value);
        }
        data.push_back(row);
    }

    int total = data.size() / 2;
    int start = std::min(skip, total);
    int count = (nquirks < 0) ? (total - start) : std::min(nquirks, total - start);
    int end = start + count;


    // after map load / bounds / KD init are done, before stepping:


    // Loop over quirks in file
    for (int h = start; h < end; ++h)
    {
        auto start = std::chrono::high_resolution_clock::now();
        // Set initial conditions
       

        int q1 = 1, q2 = -1;
        double mq = round(data[2 * h][6]); // Quirk mass in GeV
        double direc = (data[2 * h][4] + data[2 * h + 1][4] > 0) ? 1 : -1;

        std::vector<double> p1 = {direc * data[2 * h][2], direc * data[2 * h][3], direc * data[2 * h][4]};
        std::vector<double> p2 = {direc * data[2 * h + 1][2], direc * data[2 * h + 1][3], direc * data[2 * h + 1][4]};

        double E1 = sqrt(mq * mq + p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
        double E2 = sqrt(mq * mq + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);



        // Compute survival probabilities before simulation begins
        double epsilons[3] = {0.07, 0.10, 0.13};
        double epsilon1 = 1e-15;          
        double l2 = 474.64 + 5.78;    // preshower1  distance from ip
        std::string model = "f31";      // fermionic colorless model

        double decay_dist[3];
        for (int i = 0; i < 3; ++i) {
            decay_dist[i] = DecayDistance(epsilons[i],epsilon1,Lambda,mq,model,p1[0], p1[1], p1[2],p2[0], p2[1], p2[2]);
        }
        
        double decay_dist_stddev = DecayStandardDeviation(Lambda,mq,p1[0], p1[1], p1[2], E1,p2[0], p2[1], p2[2], E2);
        
        
        double survival_prob1 = LifetimeSurvivalProbGauss(decay_dist[0],decay_dist_stddev,l2);
        double survival_prob2 = LifetimeSurvivalProbGauss(decay_dist[1],decay_dist_stddev,l2);
        double survival_prob3 = LifetimeSurvivalProbGauss(decay_dist[2],decay_dist_stddev,l2);

        // Per-event transverse smearing widths at tracker plane (micrometers), for the three epsilon points
        double sigmaR[3];
        for (int i = 0; i < 3; ++i) {
            sigmaR[i] = SigmaR_um(epsilons[i], l2, Lambda, mq,
                                  p1[0], p1[1], p1[2],
                                  p2[0], p2[1], p2[2]);
        }


        std::vector<double> v1 = {p1[0] / E1, p1[1] / E1, p1[2] / E1};
        std::vector<double> v2 = {p2[0] / E2, p2[1] / E2, p2[2] / E2};

        double v10 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        double v20 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);

        std::vector<double> Beta = {(p1[0] + p2[0]) / (E1 + E2), (p1[1] + p2[1]) / (E1 + E2), (p1[2] + p2[2]) / (E1 + E2)};
        std::vector<long double> p1_long(3);
        std::vector<long double> p2_long(3);
        long double E1_long = static_cast<double long>(E1);
        long double E2_long = static_cast<double long>(E2);

        for (int i = 0; i < p1.size(); ++i)
        {
            p1_long[i] = static_cast<long double>(p1[i]);
            p2_long[i] = static_cast<long double>(p2[i]);
        }

        std::vector<long double> Beta_long = {(p1_long[0] + p2_long[0]) / (E1_long + E2_long), (p1_long[1] + p2_long[1]) / (E1_long + E2_long), (p1_long[2] + p2_long[2]) / (E1_long + E2_long) };

        long double mq_long = static_cast<long double>(mq);
        long double Lambda_long = static_cast<long double>(Lambda);

        long double t1q_long = 658 * ((2 * mq_long) / (Lambda_long * Lambda_long)) * sqrt(pow((E1_long + E2_long) / (2 * mq_long), 2) - 1 / (1 - (Beta_long[0] * Beta_long[0] + Beta_long[1] * Beta_long[1] + Beta_long[2] * Beta_long[2])));

        double t1q = static_cast<double>(t1q_long);

       
        double beta_mag = std::sqrt(Beta[0]*Beta[0] + Beta[1]*Beta[1] + Beta[2]*Beta[2]);
        const double c = 2.99792458e8;  
        double distance_period = beta_mag * c * t1q * 1e-9;  // meters
        std::cout<<"closest approach interval: " << distance_period << " m" << std::endl;
        // Decide whether we're using the "within half-oscillation" shortcut between front and back.
        double half_osc_length = distance_period * 1e6; // (meters -> micrometers) for closest-approach interval
        const double eps = 1e-9;
        bool within_half = (back - front) <= (half_osc_length * (1.0 + eps));

        // ------------------------------------------------------------------
        // Phase mapping / initialization
        // We first compute the original (pre-loss) phase-mapped start position
        // (r1=r2) exactly as before, then (if within_half) we apply ionization
        // loss to momenta and recompute ONLY the start times t1/t2 using the
        // post loss Beta_z and t1q, while keeping r1/r2 fixed 
        // ------------------------------------------------------------------

        // Original pre-loss phase map to the start plane (near z=front)
        int nsf0 = floor(front / (3e5 * t1q * Beta[2]));
        int ns0 = (nsf0 % 2 == 0) ? nsf0 : nsf0 - 1;
        double t1_pre = ns0 * t1q;
        double t2_pre = t1_pre;
        std::vector<double> r1_pre = {3e5 * ns0 * t1q * Beta[0], 3e5 * ns0 * t1q * Beta[1], 3e5 * ns0 * t1q * Beta[2]};
        std::vector<double> r2_pre = r1_pre;
        const double z_start = r1_pre[2];

        // If we start the simulation at `front` near `back`, apply a fast ionization-loss correction for the skipped 0 -> front segment.
        if (within_half && front > 0.0) {
            const int mq_int = static_cast<int>(std::llround(mq));
            auto recompute_kin = [&]() {
                v1 = DivideVector(p1, E1);
                v2 = DivideVector(p2, E2);
                v10 = std::sqrt(DotProduct(v1, v1));
                v20 = std::sqrt(DotProduct(v2, v2));
                Beta = {(p1[0] + p2[0]) / (E1 + E2), (p1[1] + p2[1]) / (E1 + E2), (p1[2] + p2[2]) / (E1 + E2)};

                for (int i = 0; i < 3; ++i) {
                    p1_long[i] = static_cast<long double>(p1[i]);
                    p2_long[i] = static_cast<long double>(p2[i]);
                }
                E1_long = static_cast<long double>(E1);
                E2_long = static_cast<long double>(E2);
                Beta_long = {(p1_long[0] + p2_long[0]) / (E1_long + E2_long),
                             (p1_long[1] + p2_long[1]) / (E1_long + E2_long),
                             (p1_long[2] + p2_long[2]) / (E1_long + E2_long)};
                mq_long = static_cast<long double>(mq);
                Lambda_long = static_cast<long double>(Lambda);
                t1q_long = 658 * ((2 * mq_long) / (Lambda_long * Lambda_long)) *
                           sqrt(pow((E1_long + E2_long) / (2 * mq_long), 2) -
                                1 / (1 - (Beta_long[0] * Beta_long[0] + Beta_long[1] * Beta_long[1] + Beta_long[2] * Beta_long[2])));
                t1q = static_cast<double>(t1q_long);
                beta_mag = std::sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2]);
                distance_period = beta_mag * c * t1q * 1e-9;  // meters
                half_osc_length = distance_period * 1e6;
                within_half = (back - front) <= (half_osc_length * (1.0 + eps));
            };

            // Apply ionization loss to momenta for the skipped 0 -> front segment.
            // (We keep r1/r2 fixed to the pre-loss placement computed above.)
            double unused_tof = 0.0;
            const bool ok = RangeTables::apply_range_table_loss(
                mq_int, mq, Lambda, front, distance_period, Beta, beta_cutOff, p1, E1, p2, E2, unused_tof
            );
            if (!ok) continue;

            // recompute post-loss Beta and t1q for the TOF/phase correction
            recompute_kin();

            // recompute the phase mapped time using the fixed start z (from the original placement).
         
            int nsf_new = floor(z_start / (3e5 * t1q * Beta[2]));
            int ns_new = (nsf_new % 2 == 0) ? nsf_new : nsf_new - 1;
            t1_pre = ns_new * t1q;
            t2_pre = t1_pre;

            
        }

        double dt = std::min(0.03, t1q / divStep);

        // Final start state:
        // - positions are the original pre-loss placement (r1_pre/r2_pre)
        // - times are corrected (t1_pre/t2_pre), using post-loss Beta_z and t1q but fixed z_start
        double t1 = t1_pre;
        double t2 = t2_pre;
        std::vector<double> r1 = r1_pre;
        std::vector<double> r2 = r2_pre;

        
        int stepcount = 0;
        int n = 1;
        double lastSaveTime = 0;
        double saveInterval = .1; // nano seconds
        double dx1pre = 0.0, dx2pre = 0.0;

        double momentum_diff_max = 0.0;  // Maximum momentum difference in COM frame
        std::deque<double> momentum_diff_window;  // Store momentum differences in COM frame
        struct Snapshot
        {
            std::vector<double> r1, r2, p1, p2;
            double t1, t2;
            int trajectory_lines;  // Number of trajectory lines written up to this point
        };
        std::deque<Snapshot> state_window;
     
        // half_osc_length / within_half already defined above (and recomputed if we applied 0->front energy loss)
        
        // Calculate the z-position of the second-to-last meeting point before back
        // Meeting points occur every osc_length, so second-to-last is at back - half_osc_length
        double window_factor = 1.5;
        double second_to_last_meeting_z = back - window_factor * half_osc_length;
        bool trackingMinimum = false;  // Flag to start tracking once past second-to-last meeting (keep name for consistency)
        int stepsSinceTrackingStart = 0;
        int deque_size = static_cast<int>(divStep * window_factor);  // Deque size based on divStep
        int trajectory_steps_written = 0;  // Track how many trajectory points we've written
        bool shouldSave = false;  // Only save if we hit back or within_half, not if too slow/transverse


        // main step loop
        while (!((sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2]) < beta_cutOff) ||
                 (sqrt(((r1[0] + r2[0]) / 2) * ((r1[0] + r2[0]) / 2) + ((r1[1] + r2[1]) / 2) * ((r1[1] + r2[1]) / 2)) > 1.5e6)))
        { // if quirks transverse cm goes a 2.0m off beamline, or too slow cancel event

            int loct1 = Loct(r1[0], r1[1], r1[2]);
            int loct2 = Loct(r2[0], r2[1], r2[2]);

            // int layer1i = Layer(r1[0], r1[1], r1[2]);
            // int layer2i = Layer(r2[0], r2[1], r2[2]);
            stepcount++;

            // Recalculate energies and velocities
            E1 = sqrt(mq * mq + p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
            v1 = DivideVector(p1, E1);
            v10 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

            E2 = sqrt(mq * mq + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);
            v2 = DivideVector(p2, E2);
            v20 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);

            Beta = {(p1[0] + p2[0]) / (E1 + E2), (p1[1] + p2[1]) / (E1 + E2), (p1[2] + p2[2]) / (E1 + E2)};

            // Direction of the string at quirk 1
            std::vector<double> s1(3), s2(3);
            if (n == 1)
            {
                s1 = Normalize(SubtractVectors(v1, v2));
            }
            else
            {
                s1 = Normalize(AddVectors(MultiplyVector(SubtractVectors(SubtractVectors(r1, r2), MultiplyVector(Beta, 300000 * (t1 - t2))), (1 - DotProduct(Beta, Beta))), MultiplyVector(SubtractVectors(Beta, v1), DotProduct(SubtractVectors(SubtractVectors(r1, r2), MultiplyVector(Beta, 300000 * (t1 - t2))), Beta))));
            }

            double vp1 = DotProduct(v1, s1);
            std::vector<double> vc1 = SubtractVectors(v1, MultiplyVector(s1, vp1));
            double vc10 = sqrt(DotProduct(vc1, vc1));

            // Direction of the string at quirk 2
            if (n == 1)
            {
                s2 = Normalize(SubtractVectors(v2, v1));
            }
            else
            {
                s2 = Normalize(AddVectors(MultiplyVector(SubtractVectors(SubtractVectors(r2, r1), MultiplyVector(Beta, 300000 * (t2 - t1))), (1 - DotProduct(Beta, Beta))), MultiplyVector(SubtractVectors(Beta, v2), DotProduct(SubtractVectors(SubtractVectors(r2, r1), MultiplyVector(Beta, 300000 * (t2 - t1))), Beta))));
            }

            double vp2 = DotProduct(v2, s2);
            std::vector<double> vc2 = SubtractVectors(v2, MultiplyVector(s2, vp2));
            double vc20 = sqrt(DotProduct(vc2, vc2));

            // Estimation of the travel distance (in cm) using the average de/dx value

            if (loct1 > 0 || loct2 > 0)
            {

                std::vector<double> F1pre = CalculateForces(mq, Lambda, v1, vc1, s1, vc10, vp1, loct1, q1, r1);
                std::vector<double> F2pre = CalculateForces(mq, Lambda, v2, vc2, s2, vc20, vp2, loct2, q2, r2);

                double ct1pre = CalculateCt(v1, Beta, r1, r2, F1pre, E1, E2);
                double ct2pre = CalculateCt(v2, Beta, r2, r1, F2pre, E1, E2);

                double dt1pre, dt2pre;
                if (abs(ct1pre) < abs(ct2pre))
                {
                    dt1pre = dt;
                    dt2pre = dt1pre * ct1pre / ct2pre;
                }
                else
                {
                    dt2pre = dt;
                    dt1pre = dt2pre * ct2pre / ct1pre;
                }

                std::vector<double> p1pre = AddVectors(p1, MultiplyVector(F1pre, dt1pre));
                std::vector<double> p2pre = AddVectors(p2, MultiplyVector(F2pre, dt2pre));

                dx1pre = CalculateDistance(v1, p1pre, mq, dt1pre);
                dx2pre = CalculateDistance(v2, p2pre, mq, dt2pre);
            }

            // Recalculate forces using normally distributed de/dx
            std::vector<double> F1 = CalculateForcesWithGaus(mq, Lambda, v1, vc1, s1, vc10, vp1, loct1, q1, r1, dx1pre, gen);
            std::vector<double> F2 = CalculateForcesWithGaus(mq, Lambda, v2, vc2, s2, vc20, vp2, loct2, q2, r2, dx2pre, gen);

            double ct1 = CalculateCt(v1, Beta, r1, r2, F1, E1, E2);
            double ct2 = CalculateCt(v2, Beta, r2, r1, F2, E1, E2);

            double dt1, dt2;

            if (abs(ct1) < abs(ct2))
            {
                dt1 = dt;
                dt2 = dt1 * ct1 / ct2;
            }
            else
            {
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

 
            if (!std::isfinite(r1[0]) || !std::isfinite(r1[1]) || !std::isfinite(r1[2]) ||
                !std::isfinite(r2[0]) || !std::isfinite(r2[1]) || !std::isfinite(r2[2])) {
                std::cerr << "NaN detected in positions at event " << h << ", skipping.\n";
                break; // skip this event
            }

        

            // Determine the detector scintillators, currently not used
            // int layer1f = Layer(r1[0], r1[1], r1[2]);
            // int layer2f = Layer(r2[0], r2[1], r2[2]);

            // write trajectory info to file
            if (traj && outputFileTrajectory && outputFileTrajectory->is_open())
            {
                if (t1 - lastSaveTime >= saveInterval)
                {
                    
                    double pair_px = p1[0] + p2[0];
                    double pair_py = p1[1] + p2[1];
                    double pair_pz = p1[2] + p2[2];
                    double pair_p_mag = sqrt(pair_px * pair_px + pair_py * pair_py + pair_pz * pair_pz);
                    double pair_speed = pair_p_mag / (E1 + E2);
                    
                    *outputFileTrajectory << std::setprecision(16) << t1 << " " << r1[0] << " " << r1[1] << " " << r1[2] << " " << pair_speed << "\n";
                    *outputFileTrajectory << std::setprecision(16) << t1 << " " << r2[0] << " " << r2[1] << " " << r2[2] << " " << pair_speed << "\n";
                    lastSaveTime = t1;
                    trajectory_steps_written += 2;  
                }
            }

            if (stepcount % 100000 == 0)
            {
                std::cout << "r1: " << r1[0] << " " << r1[1] << " " << r1[2] << std::endl;
                std::cout << "r2: " << r2[0] << " " << r2[1] << " " << r2[2] << std::endl;
                std::cout << "beta: " << sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2]) << std::endl;
            }


            // Calculate momentum difference in COM frame
            // Boost momenta to COM frame and calculate |p1_COM - p2_COM|
            std::vector<double> p1_COM = BoostToCOM(p1, E1, Beta);
            std::vector<double> p2_COM = BoostToCOM(p2, E2, Beta);
            std::vector<double> p_diff = SubtractVectors(p1_COM, p2_COM);
            double momentum_diff = sqrt(DotProduct(p_diff, p_diff));  // |p1_COM - p2_COM|

            // Check if we've passed the second-to-last meeting point and start tracking (using OR)
            if (!trackingMinimum && (r1[2] >= second_to_last_meeting_z || r2[2] >= second_to_last_meeting_z))
            {
                trackingMinimum = true;
                momentum_diff_max = momentum_diff;
            }

            if (trackingMinimum)
            {
                momentum_diff_window.push_back(momentum_diff);
                state_window.push_back({r1,r2,p1,p2,t1,t2,trajectory_steps_written});

                if (momentum_diff_window.size() > deque_size) {
                    momentum_diff_window.pop_front();
                    state_window.pop_front();
                }
            }
            
            // Update maximum momentum difference while tracking
            if (trackingMinimum)
            {
                if (momentum_diff > momentum_diff_max)
                {
                    momentum_diff_max = momentum_diff;
                }
            }

            // Only break if within_half (separate pathway) - output done after loop
            if (within_half)
            {
                shouldSave = true;
                break;
            }
            
            // Break if we hit back (not within_half) - will search deque after loop
            if (!within_half && (r1[2] >= back || r2[2] >= back))
            {
                shouldSave = true;
                break;
            }

            n++;
        }
        
        // After loop: if we hit back (not within_half), search deque for maximum momentum difference
        Snapshot maxSnapshot;
        if (!within_half && (r1[2] >= back || r2[2] >= back) && trackingMinimum && momentum_diff_window.size() > 0)
        {
            // Find maximum in deque
            int maxIdx = 0;
            double maxMomentumDiff = momentum_diff_window[0];
            for (size_t i = 1; i < momentum_diff_window.size(); ++i) {
                if (momentum_diff_window[i] > maxMomentumDiff) {
                    maxMomentumDiff = momentum_diff_window[i];
                    maxIdx = i;
                }
            }
            maxSnapshot = state_window[maxIdx];
            momentum_diff_max = maxMomentumDiff;
            
            // Restore to maximum momentum difference position (closest approach)
            r1 = maxSnapshot.r1;
            r2 = maxSnapshot.r2;
            p1 = maxSnapshot.p1;
            p2 = maxSnapshot.p2;
            t1 = maxSnapshot.t1;
            t2 = maxSnapshot.t2;
            
            // Truncate trajectory file to the point where maximum occurred
            if (traj && outputFileTrajectory) {
                if (outputFileTrajectory->is_open()) {
                    outputFileTrajectory->close();
                }
                
                std::string trajFileName = stem + "_" + lambdaStr + "eV_trajectory.txt";
                std::ifstream trajIn(trajFileName);
                std::vector<std::string> lines;
                std::string line;
                while (std::getline(trajIn, line)) {
                    lines.push_back(line);
                }
                trajIn.close();
                
                // Keep only the first trajectory_lines from the maximum snapshot
                int lines_to_keep = maxSnapshot.trajectory_lines;
                if (lines_to_keep > 0 && lines_to_keep <= static_cast<int>(lines.size())) {
                    lines.erase(lines.begin() + lines_to_keep, lines.end());
                }
                
                // Write back truncated file
                std::ofstream trajOut(trajFileName);
                for (const auto& l : lines) {
                    trajOut << l << "\n";
                }
                trajOut.close();
            }
        }


        double t_star = (t1 > t2) ? t1 : t2;
        SyncQuirksToSameTime(t_star, r1, p1, t1, q1, r2, p2, t2, q2, mq, Lambda);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        // Calculate distance at closest approach point
        double dist_at_closest = sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]));
        std::cout << "Closest approach found at z = " << r1[2] << " (r1), " << r2[2] << " (r2) with distance: " << dist_at_closest << " and momentum diff (COM): " << momentum_diff_max << std::endl;
        std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;


        
        if (shouldSave){
        // Sync and output (for both within_half and deque cases)
    
        outputFile << std::setprecision(16) << h << " " << mq << " " << Lambda << " " << t1q << " 1 " << 1 << " " << t1 << " "
                   << r1[0] << " " << r1[1] << " " << r1[2] << " " << p1[0] << " " << p1[1] << " " << p1[2] << " " << duration.count()
                   << " " << survival_prob1 << " " << survival_prob2 << " " << survival_prob3
                   << " " << sigmaR[0] << " " << sigmaR[1] << " " << sigmaR[2] << "\n";
        outputFile << std::setprecision(16) << h << " " << mq << " " << Lambda << " " << t1q << " 2 " << 1 << " " << t2 << " "
                   << r2[0] << " " << r2[1] << " " << r2[2] << " " << p2[0] << " " << p2[1] << " " << p2[2] << " " << duration.count()
                   << " " << survival_prob1 << " " << survival_prob2 << " " << survival_prob3
                   << " " << sigmaR[0] << " " << sigmaR[1] << " " << sigmaR[2] << "\n";
        }
        
        double final_beta = sqrt(Beta[0] * Beta[0] + Beta[1] * Beta[1] + Beta[2] * Beta[2]);
        double final_trans = sqrt(((r1[0] + r2[0]) / 2) * ((r1[0] + r2[0]) / 2) + ((r1[1] + r2[1]) / 2) * ((r1[1] + r2[1]) / 2));
        std::cout << "DEBUG: Loop exited. beta=" << final_beta << " (cutoff=" << beta_cutOff << "), trans_dist=" << final_trans << " (limit=1.5e6)" << std::endl;
        std::cout << h << std::endl;
      
    }

    inputFile.close();
    outputFile.close();

    if (traj && outputFileTrajectory )
    {
        outputFileTrajectory->close();
    }

    return 0;
}
