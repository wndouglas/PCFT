#include "tools/Timer.hpp"
#include "catch.hpp"
#include <thread>

using namespace PCFT::tools;
using namespace std;

TEST_CASE( "Test Timer accuracy", "[timer]" )
{
    // Arrange
    Timer timer;

    // Act
    using namespace std::chrono_literals;
    auto pauseTime = 50ms;
    Timer::milliseconds pauseTimeMs = pauseTime.count();

    // Best of three calculations
    const size_t N = 3;
    Timer::milliseconds ave = 0.0;
    for(size_t n = 0; n < N; n++)
    {
        timer.start();
        this_thread::sleep_for(pauseTime);
        timer.stop();
        ave += timer.duration()/N;
    }

    // Assert
    REQUIRE(fabs(ave - pauseTimeMs)/pauseTimeMs < 0.1); // Relative error within 10%
}