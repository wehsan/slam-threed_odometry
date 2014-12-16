#include <boost/test/unit_test.hpp>
#include <3d_odometry/Dummy.hpp>

using namespace 3d_odometry;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    3d_odometry::DummyClass dummy;
    dummy.welcome();
}
