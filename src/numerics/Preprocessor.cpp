#include "numerics/Preprocessor.hpp"

using namespace PCFT::numerics;
using namespace std;

Preprocessor::Preprocessor(std::unique_ptr<IFourierTransformer> transformer) : mTransformer(std::move(transformer)) { }

void Preprocessor::execute(const vector<double>& inputVector, vector<double>& outputVector)
{
    //mTransformer->fft()
}