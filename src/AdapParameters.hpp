/* This library implements the Adaptive Parameters Identification method
 *
 * Smallwood, D. A., & Whitcomb, L. L. (2003). Adaptive identification of
 * dynamically positioned underwater robotic vehicles.
 *
 * It uses the simplified motion model, where the all 6 degree of freedom
 * are decoupled, the auv is symmetrical in its 3 planes, the auv performs at
 * low speed, with modest attitude changes and with single-degree-of-freedom
 * motion (auv moving at one direction at a time).
 *
 * That means, the equation for each degree of freedom will be of the form:
 * taui(t)= mi*vdoti + dQi*vi(t)*|vi(t)| + dLi*vi(t) + bi,
 * i=1...6, mi,dQi,dLi>0
 */

#ifndef ADAP_PARAMETERS_HPP_
#define ADAP_PARAMETERS_HPP_

#include "DataTypes.hpp"
#include <iostream>
#include <math.h>
#include <vector>

namespace adap_parameters_estimator{

class AdapParameters {

public:

    AdapParameters(base::Vector4d gain_lambda, double gain_a, double step = 0.01, double maximum_inertia = 10000, double initial_inertia = 100);
    ~AdapParameters();

    /** Set gains of adaptive identifier
     *
     * Expected valeus are gain_a<0, gain_lambda>0
     * @param gain_lambda, for the parameters
     * @param gain_a, for the velocity error
     */
    void setGains(base::Vector4d gain_lambda, double gain_a);

    /** Set step
     *
     * @param step
     */
    void setStep(double step);

    /** Convert the parameters used in the adaptive law for those used in the dynamic model
     *
     * @param phi_parameters, parameters used in adaptive law
     * @return parameters of dynamic model: [m, dq, dl, b]
     */
    base::Vector4d convertParameters(const base::Vector4d &phi_parameters);

    /** Reset internal parameters
     */
    void resetStates(void);

    /** Set initial parameters
     *
     *  @param parameters, [inertia, quad_damp, lin_damp, buoy]
     */
    void setParameters(const base::Vector4d &parameters);

    /** Estimate parameters based on adaptive law
     *
     * @param effort, force or torque applied in one DOF
     * @param velocity measured in one DOF
     * @return vector with estimated parameters
     */
    base::Vector4d estimateParameters(double effort, double velocity);

    /** Get delta Velocity
     *
     * @return _delta_velocity
     */
    double getDeltaVelocity(void);


private:
    base::Vector4d _estimated_phi;
    double _estimated_velocity;
    double _delta_velocity;
    double _gain_a;
    base::Vector4d _gain_lambda;
    double _step;
    double _maximum_inertia;
    double _minimum_inertia;
};
}

#endif /* ADAP_PARAMETERS_HPP_ */
