#include "AdapParameters.hpp"

namespace adap_parameters_estimator{

AdapParameters::AdapParameters(base::Vector4d gain_lambda, double gain_a, double step, double maximum_inertia, double initial_inertia)
{
    resetStates();
    setGains(gain_lambda, gain_a);
    setStep(step);
    if(maximum_inertia <= 0)
        throw std::runtime_error("AdapParameters: maximum_inertia must have a positive value");
    _maximum_inertia = maximum_inertia;
    _minimum_inertia = initial_inertia;
    setParameters(base::Vector4d(initial_inertia,0,0,0));
}

AdapParameters::~AdapParameters()
{
}

void AdapParameters::setGains(base::Vector4d gain_lambda, double gain_a)
{
    if(gain_a >= 0)
        throw std::runtime_error("AdapParameters: checkGains: The gain_a must have a negative value");
    for(uint i=0; i<4; i++)
    {
        if (gain_lambda[i] <= 0)
	       throw std::runtime_error("AdapParameters: checkGains: The gain_lambda must have a positive value");
    }
    _gain_a = gain_a;
    _gain_lambda = gain_lambda;
}

void AdapParameters::setStep(double step)
{
    if(step <= 0)
        throw std::runtime_error("AdapParameters: setStep: Step must have a positive value");
    _step = step;
}

base::Vector4d AdapParameters::convertParameters(const base::Vector4d &phi_parameters)
{
    base::Vector4d parameters;
    double alpha = phi_parameters[0];
    if (alpha < 1/_maximum_inertia)
        alpha = 1/_maximum_inertia;

    parameters[0] = 1/alpha;
    parameters[1] = -phi_parameters[1]/alpha;
    parameters[2] = -phi_parameters[2]/alpha;
    parameters[3] = -phi_parameters[3]/alpha;
    return parameters;
}

void AdapParameters::resetStates(void)
{
    _estimated_phi = base::Vector4d(1/_minimum_inertia,0,0,0);
    _estimated_velocity = 0;
    _delta_velocity = 0;
    _gain_a = -1;
    _gain_lambda = base::Vector4d(1,1,1,1);
    _step = 0.01;
}

void AdapParameters::setParameters(const base::Vector4d &parameters)
{
    if(parameters[0] <= 0)
        throw std::runtime_error("AdapParameters: setParameters: inertia must have a positive value");
    _estimated_phi[0] = 1/parameters[0];
    _estimated_phi[1] = -parameters[1]/parameters[0];
    _estimated_phi[2] = -parameters[2]/parameters[0];
    _estimated_phi[3] = -parameters[3]/parameters[0];
}

base::Vector4d AdapParameters::estimateParameters(double effort, double velocity)
{
    // Error of velocity
    _delta_velocity = _estimated_velocity - velocity;
    // Estimated states
    base::Vector4d states(effort, (_estimated_velocity*fabs(_estimated_velocity)), _estimated_velocity, 1);
    // Adaptive estimator
    double acceleration = (_gain_a * _delta_velocity) + (_estimated_phi.transpose() * states);
    // Update law
    base::Vector4d phi_dot = -1*_gain_lambda.asDiagonal() * _delta_velocity * states;

    // Euler integrator of velocity
    _estimated_velocity += acceleration * _step;
    // Euler integrator of parameters
    _estimated_phi += phi_dot * _step;

    // Transform into conventional parameters (inertia, damping, buoyancy)
    return convertParameters(_estimated_phi);
}

double AdapParameters::getDeltaVelocity()
{
    return _delta_velocity;
}

};
