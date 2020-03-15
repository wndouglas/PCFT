// PCFT.cpp : Defines the entry point for the application.
//

#include "tools/Timer.hpp"
#include "numerics/FTFactory.hpp"
#include "numerics/Preprocessor.hpp"
#include "DomainParameters.hpp"
#include "PCFTExecutor.hpp"
#include <complex>
#include <cmath>
#include <thread>
#include <iostream>
#include <fstream>
#include <string>

using namespace PCFT;
using namespace std;

void do_transform(string file_name, int N, int M, double r, double sigma, double T, double Smin, double Smax, double eps1, double eps2);
void output_to_file(string file_name, const OutputPage& op, const DomainParameters& inputPage, tools::Timer::milliseconds dur, const double r, const double sigma, const double T);

int main(int argc, const char * argv[])
{
    std::string file_name;
    int N;
    int M;
    double r;
    double sigma;
    double T;
    double Smin;
    double Smax;
    double epsilon1;
    double epsilon2;
    try
    {
        file_name = argv[1];
        N = std::stoi(argv[2]);
        M = std::stoi(argv[3]);
        r = std::stod(argv[4]);
        sigma = std::stod(argv[5]);
        T = std::stod(argv[6]);
        Smin = std::stod(argv[7]);
        Smax = std::stod(argv[8]);
        epsilon1 = std::stod(argv[9]);
        epsilon2 = std::stod(argv[10]);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cerr << "Given file name: " << file_name << '\n';
        std::cerr << "N = " << N << '\n';
        std::cerr << "M = " << M << '\n';
        std::cerr << "r = " << r << '\n';
        std::cerr << "sigma = " << sigma << '\n';
        std::cerr << "T = " << T << '\n';
        std::cerr << "Smin = " << Smin << '\n';
        std::cerr << "Smax = " << Smax << '\n';
        std::cerr << "epsilon1 = " << epsilon1 << '\n';
        std::cerr << "epsilon2 = " << epsilon2 << '\n';
    }

    do_transform(file_name, N, M, r, sigma, T, Smin, Smax, epsilon1, epsilon2);

	return 0;
}

void do_transform(string file_name, int N, int M, double r, double sigma, double T, double Smin, double Smax, double eps1, double eps2)
{
    using namespace numerics;
    using namespace tools;

    double xMin = std::log(Smin);
    double xMax = std::log(Smax);

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        M,
        N,

        xMin,
        xMax,
        T,

        eps1,
        eps2
    };
    
    PCFTExecutor executor(
        FTFactory::instance(N*2),
        std::make_unique<Preprocessor>(GFunction(r, sigma, DomainParameters::getDTau(pPackage.T, pPackage.M))),
        pPackage);

    OutputPage outPage;
    Timer t;
    t.start();
    executor.execute(outPage);
    t.stop();
    Timer::milliseconds dur = t.duration();

    output_to_file(file_name, outPage, pPackage, dur, r, sigma, T);
}

void output_to_file(string file_name,
    const OutputPage& op,
    const DomainParameters& pPackage, 
    tools::Timer::milliseconds dur,
    const double r, 
    const double sigma, 
    const double T)
{
    ofstream myfile;
    myfile.open (file_name);

    myfile << "S, V(S, 0), M, N, x_min, x_max, r, sigma, T, CalcTime" << std::endl;
    for (int i = 0; i < pPackage.N; i++)
    {
        myfile << op.state[i] << ", " << op.valueFunction[0][i];
        if(i == 0)
            myfile << ", " << pPackage.M << ", " << pPackage.N << ", " << pPackage.xMin << ", " << pPackage.xMax << ", " << r << ", " << sigma << ", " << T << ", " << dur;
        else
            myfile << ", 0, 0, 0, 0, 0, 0, 0, 0";

        myfile << std::endl;
    }

    myfile.close();
}

