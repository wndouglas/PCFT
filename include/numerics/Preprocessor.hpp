#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <vector>
#include "numerics/IFourierTransformer.hpp"
#include "ParameterPackage.hpp"

namespace PCFT
{
	namespace numerics
	{
		class Preprocessor
		{
		public:
			Preprocessor(
				std::unique_ptr<IFourierTransformer> transformer,
				ParameterPackage pPackage);

			void execute(const IFourierTransformer::RealVec& inputVector, IFourierTransformer::RealVec& outputVector) const;

		private:
			std::unique_ptr<IFourierTransformer> mTransformer;
			const int mN;
			const double mDx;
			const double mDtau;
			const double mEpsilon1;
			const double mEpsilon2;

			void shiftedFft(IFourierTransformer::ComplexVec& inputVector, IFourierTransformer::RealVec& outputVector) const;
			void shiftedIfft(const IFourierTransformer::RealVec& inputVector, IFourierTransformer::ComplexVec& outputVector) const;
		};
	}
}

#endif // !PREPROCESSOR_HPP