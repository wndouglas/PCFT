#include "tools/Timer.hpp"

using namespace PCFT::tools;
using namespace std;

Timer::Timer() : mClockRunning(false), mHasClockStopped(false)
{
    mStartTime = chrono::high_resolution_clock::now();
    mFinishTime = mStartTime;
}

void Timer::start()
{
    mStartTime = chrono::high_resolution_clock::now();
    mClockRunning = true;
}

void Timer::stop()
{
    mFinishTime = chrono::high_resolution_clock::now();
    mClockRunning = false;
    mHasClockStopped = true;
}

Timer::milliseconds Timer::duration()
{
    if(!mHasClockStopped)
    {
        return chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - mStartTime).count();
    }

    double timeTaken;
    if(!mClockRunning)
    {
        timeTaken = chrono::duration_cast<chrono::milliseconds>(mFinishTime - mStartTime).count();
        mStartTime = mFinishTime;
    }
    else
    {
        auto currTime = chrono::high_resolution_clock::now();
        timeTaken = chrono::duration_cast<chrono::milliseconds>(currTime - mStartTime).count();
    }
    return timeTaken;
}