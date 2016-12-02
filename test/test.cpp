#define BOOST_TEST_MODULE ADAP_PARAMETERS_ESTIMATOR
#include <adap_parameters_estimator/AdapParameters.hpp>
#include <boost/test/included/unit_test.hpp>

using namespace adap_parameters_estimator;

BOOST_AUTO_TEST_SUITE (CONSTRUCTOR)

BOOST_AUTO_TEST_CASE( normal_constructor)
{
    base::Vector4d gain_lambda(1,1,1,1);
    double gain_a = -1;
    double step = 0.01;
	BOOST_REQUIRE_NO_THROW( AdapParameters adaptive(gain_lambda, gain_a, step));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE (METHODS)

BOOST_AUTO_TEST_CASE( get_parameters)
{
    base::Vector4d gain_lambda(1,1,1,1);
    double gain_a = -1;
    double step = 0.01;
    AdapParameters adaptive(gain_lambda, gain_a, step);
    gain_lambda = base::Vector4d(2,2,2,2);
    gain_a = -2;
    BOOST_REQUIRE_NO_THROW( adaptive.setGains(gain_lambda, gain_a));
    BOOST_REQUIRE_NO_THROW( adaptive.setStep(0.2));
    BOOST_REQUIRE_NO_THROW( adaptive.resetStates());
    BOOST_REQUIRE_EQUAL( adaptive.getDeltaVelocity(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
