#include <boost/test/unit_test.hpp>
#include <adap_parameters_estimator/Dummy.hpp>

using namespace adap_parameters_estimator;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    adap_parameters_estimator::DummyClass dummy;
    dummy.welcome();
}
