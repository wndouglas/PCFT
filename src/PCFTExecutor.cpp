#include "PCFTExecutor.hpp"

using namespace PCFT;

void PCFTExecutor::execute(OutputPage& outputPage)
{
    DomainParameters mExtendedParams
    {
        mParams.M,
        mParams.N*2,

        mParams.xMin - 0.5*(mParams.xMax - mParams.xMin),
        mParams.xMax + 0.5*(mParams.xMax - mParams.xMin),
        mParams.T,

        mParams.epsilon1,
        mParams.epsilon2
    };

    numerics::IFourierTransformer::ComplexVec GOut = mPreprocessor->execute(mTransformer.get(), mExtendedParams);

    // For now we try a European call
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;

    double xMin = mParams.xMin;
    double xMax = mParams.xMax;
    double N = mParams.N;
    double M = mParams.M;
    double NExtended = mExtendedParams.N;

    double dx = DomainParameters::getDx(xMax, xMin, N);
    std::vector<double> SGrid(N);
    std::vector<double> SGridExtended(NExtended);
    std::vector<double> vGrid(N);
    std::vector<double> vGridExtended(NExtended);

    double x = xMin;
    for(int i = 0; i < N; i++)
    {
        double S = std::exp(x);
        SGrid[i] = S;
        vGrid[i] = std::max(0.0, S - K);
        x += dx;
    }

    x = mExtendedParams.xMin;
    for(int i = 0; i < NExtended; i++)
    {
        double S = std::exp(x);
        vGridExtended[i] = std::max(0.0, S - K);
        SGridExtended[i] = S;
        x += dx;
    }
    
    numerics::IFourierTransformer::ComplexVec vTempGrid(NExtended);
    for(int i = 0; i < NExtended; i++)
    {
        vTempGrid[i].real(vGridExtended[i]);
    }

    for(int m = 0; m < M; m++)
    {
        mTransformer->shiftedFft(vTempGrid, vTempGrid);
        for(int i = 0; i < NExtended; i++)
        {
            vTempGrid[i] = vTempGrid[i] * GOut[i];
        }
        mTransformer->shiftedIfft(vTempGrid, vTempGrid);

        for(int i = 0; i < NExtended; i++)
        {
            if(i < N/2)
            {
                vTempGrid[i] = vTempGrid[N/2 + 1];
            }
            else if(i > 3*N/2)
            {
                vTempGrid[i] = SGridExtended[i] - K*exp(-r*T);
            }
        }
    }

    std::vector<std::vector<double> >& valueFunction = outputPage.valueFunction;
    std::vector<double>& stateVector = outputPage.state;
    valueFunction.resize(1);
    for(int i = 0; i < N; i++)
    {
        valueFunction[0].push_back(vTempGrid[i + N/2].real());
        stateVector.push_back(SGrid[i]);
    }
}