#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <vector>
#include "IFourierTransformer.hpp"

namespace PCFT
{
	namespace numerics
	{
		class Preprocessor
		{
		public:
			Preprocessor(std::unique_ptr<IFourierTransformer> transformer);
			void execute(const std::vector<double>& inputVector, std::vector<double>& outputVector);
		private:
			std::unique_ptr<IFourierTransformer> mTransformer;
		};
	}
}

#endif // !PREPROCESSOR_HPP