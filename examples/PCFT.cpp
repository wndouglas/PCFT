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

void do_transform(string file_name);
void output_to_file(string file_name, const OutputPage& op, const DomainParameters& inputPage, tools::Timer::milliseconds dur, const double r, const double sigma, const double T);

int main(int argc, const char * argv[])
{
    std::string file_name = "PCFT_out.txt";
    if (argc == 2)
        file_name = argv[1];

    do_transform(file_name);

	return 0;
}

void do_transform(string file_name)
{
    using namespace numerics;
    using namespace tools;

    int N = 100;
    int M = 10;

    // // Parameters
    // double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;

    double xMin = std::log(0.1);
    double xMax = std::log(200);

    // Test usage of preprocessor
    DomainParameters pPackage 
    {
        M,
        N,

        xMin,
        xMax,
        T,

        1e-8,
        1e-8
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
    myfile << "-------------------- Input and parameters --------------------" << std::endl;
    myfile << "M = " << pPackage.M 
            << ", N = " << pPackage.N 
            << ", x_min = " << pPackage.xMin 
            << ", x_max = " << pPackage.xMax 
            << std::endl;
    myfile << "r = " << r << ", sigma = " << sigma << ", T = " << T << std::endl;
    myfile << std::endl;

    myfile << "-------------------- Output --------------------" << std::endl;
    myfile << "Calculation time: " << dur << std::endl;
    myfile << "S, V(S, 0)" << std::endl;

    for (int i = 0; i < pPackage.N; i++)
    {
        myfile << op.state[i] << ", " << op.valueFunction[0][i] << std::endl;
    }

    myfile.close();
}

