#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <vector>
#include "numerics/IFourierTransformer.hpp"
#include "DomainParameters.hpp"
#include "GFunction.hpp"

namespace PCFT
{
	namespace numerics
	{
		class Preprocessor
		{
		public:
			Preprocessor(
				std::unique_ptr<IFourierTransformer> transformer,
				DomainParameters pPackage);

			IFourierTransformer::RealVec execute(const GFunction& GreenTransformFunction) const;

		private:
			std::unique_ptr<IFourierTransformer> mTransformer;
			const int mN;
			const int mP;
			int mLambda;
			const double mXMin;
			const double mXMax;
			const double mDx;
			const double mDtau;
			const double mEpsilon1;
			const double mEpsilon2;

			void shiftedFft(IFourierTransformer::ComplexVec& inputVector, IFourierTransformer::ComplexVec& outputVector) const;
			void shiftedIfft(const IFourierTransformer::ComplexVec& inputVector, IFourierTransformer::ComplexVec& outputVector) const;
			void calculateH(IFourierTransformer::ComplexVec& vector, const GFunction& greensFunctionTransform) const;
		};
	}
}

#endif // !PREPROCESSOR_HPP