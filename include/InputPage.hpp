#ifndef INPUT_PAGE_HPP
#define INPUT_PAGE_HPP

#include <vector>

namespace PCFT
{
    class OutputPage
    {
    public:
        OutputPage() : mUseExtendedGrid(true), mOutputPolicy(false), mOutputEntireGrid(false) {}

        OutputPage(bool useExtendedGrid, bool outputPolicy, bool outputEntireGrid) :
            mUseExtendedGrid(useExtendedGrid),
            mOutputPolicy(outputPolicy),
            mOutputEntireGrid(outputEntireGrid) {}

        bool isUsingExtendedGrid() { return mUseExtendedGrid; }
        bool isOutputtingPolicy() { return mOutputPolicy; }
        bool isOutputtingEntireGrid() { return mOutputEntireGrid; }

    private:
        bool mUseExtendedGrid;
        bool mOutputPolicy;
        bool mOutputEntireGrid;
    };
}
#endif // !INPUT_PAGE_HPP