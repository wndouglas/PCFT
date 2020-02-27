#include "numerics/Preprocessor.hpp"

using namespace PCFT::numerics;
using namespace std;

Preprocessor::Preprocessor(std::unique_ptr<IFourierTransformer> transformer,
				const int N,
				const double dx, 
				const double dTau, 
				const double epsilon1,
				const double epsilon2) :
                    mTransformer(std::move(transformer)), mN(N), 
                    mDx(dx), mDtau(dTau), mEpsilon1(epsilon1), mEpsilon2(epsilon2) { }

void Preprocessor::execute(const vector<double>& inputVector, vector<double>& outputVector)
{
    //mTransformer->fft()
}