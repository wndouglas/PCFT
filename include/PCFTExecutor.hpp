#ifndef PCFT_EXECUTOR_HPP
#define PCFT_EXECUTOR_HPP

#include "OutputPackage.hpp"
#include "numerics/IFourierTransformer.hpp"

namespace PCFT
{
    class PCFTExecutor
    {
    public:
        PCFTExecutor(
            std::unique_ptr<numerics::IFourierTransformer> transformer,
            std::unique_ptr<numerics::Preprocessor> preprocessor,
            DomainParameters pPackage) :
                mTransformer(std::move(transformer)),
                mPreprocessor(std::move(preprocessor)),
                mParams(pPackage) {}

        void execute(OutputPackage& outputPackage)
        {
            numerics::IFourierTransformer::ComplexVec GOut = mPreprocessor->execute(mTransformer, mParams);

            // For now we try a European call
            double K = 100.0;
            double T = 1.0;
            double r = 0.05;

            double xMin = mParams.xMin;
            double xMax = mParams.xMax;
            double xMinExtended = xMin - 0.5*(xMax - xMin);
            double xMaxExtended = xMax + 0.5*(xMax - xMin);
            double PExtended = xMaxExtended - xMinExtended;
            int N = mParams.N;
            int M = mParams.M;
            int NExtended = 2*N;

            double dx = DomainParameters::getDx(xMax, xMin, N);
            std::vector<double> SGrid(N);
            std::vector<double> SGridExtended(NExtended);
            std::vector<double> vGrid(N);
            std::vector<double> vGridExtended(NExtended);

            double x = mParams.xMin;
            for(int i = 0; i < N; i++)
            {
                double S = std::exp(x);
                SGrid[i] = S;
                vGrid[i] = std::max(0.0, S - K);
                x += dx;
            }

            x = xMinExtended;
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

            std::vector<std::vector<double> >& valueFunction = outputPackage.mValueFunction;
            valueFunction.resize(1);
            for(int i = 0; i < N; i++)
            {
                valueFunction[0].push_back(vTempGrid[i + N/2].real());
            }
        }


    private:
        std::unique_ptr<numerics::IFourierTransformer> mTransformer;
        std::unique_ptr<numerics::Preprocessor> mPreprocessor;
        DomainParameters mParams;
    };
}


#endif // !PCFT_EXECUTOR_HPP