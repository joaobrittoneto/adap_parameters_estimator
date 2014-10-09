/****************************************************************************/
/*  Adaptive Parameters Identification Library for AUVs	              		*/
/*                                                                         	*/
/* FILE --- adap_parameters.hpp		                           	         	*/
/*                                                                         	*/
/* PURPOSE --- Header file for a Adaptive Parameters Identification to be   */
/* used in a motion model of a underwater vehicle                          */
/*             Based on D. A. Smallwood(2003)                               */
/*                                                        	                */
/*  João Da Costa Britto Neto                                               */
/*  joao.neto@dfki.de                                                    	*/
/*  DFKI - BREMEN 2014                                                     	*/
/****************************************************************************/
/* This library implements the Adaptive Parameters Identification method
 *  based on D. A. Smallwood(2003)
 *
 * It uses the simplified motion model, where the all 6 degree of freedom
 *  are decoupled, the auv is symmetrical in its 3 planes, the auv performs at low speed, with modest attitude changes
 *  and with single-degree-of-freedom motion (auv moving at one direction at a time)
 *
 *  That means, the equation for each degree of freedom will be of the form:
 *   taui(t)= mi*vdoti + dQi*vi(t)*|vi(t)| + dLi*vi(t) + bi,  i=1...6, mi,dQi,dLi>0
 *
 */

#ifndef ADAP_PARAMETERS_HPP_
#define ADAP_PARAMETERS_HPP_


#include <eigen3/Eigen/Dense>
#include <vector>
#include "adap_dataTypes.h"
#include <iostream>
#include <queue>
#include <math.h>


namespace adap_parameters_estimator{

class AdapParameters {

public:

	AdapParameters(Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau);
	AdapParameters();
	~AdapParameters();

	void configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau);

	//Construct the vector of non-linear states
	void estimated_state (base::Vector6d tau);

	//Define deltaVelocity
	void delta_velocity (base::Vector6d _velocity);

	//Implement the adaptive identifier, which output is the v_est_dot
	double adaptive_estimator (void);

	//Update law of the parameters, which output is PHI_est_dot
	base::Vector4d update_law (void);

	//Convert the parameters used in the adaptive law for those used in the motion model
	void convetional_parameters (void);

	// Verify is the gains would lead to stability (gain_a<0, gain_lambda>0)
	void check_gains (void);

	//Euler method for establish the velocity in the online method
	void euler_velocity (void);

	//Euler method for establish the parameters in the online method
	void euler_parameters (void);

	//Construct the diagonal matrix of gain lambda for the specified dof.
	base::Matrix4d matrix_lambda (void);

	//Initial parameters values
	//void init_param (void);

	//return the conventional parameters;
	base::Vector4d get_parameters (void);

	//return the filtered conventional parameters;
	base::Vector4d get_filtered_parameters (void);

	//return deltaVelocity
	double get_delta_v (void);

	//return normDeltaVelocity
	double get_norm_delta_v (void);

	//return the non-linear vector state
	base::Vector4d get_est_states (void);

	//return the estimatedPHI
	base::Vector4d get_est_phi(void);

	//return variables used in the model
	double get_est_velocity(void);

	//return frequency of thruster input signal
	double get_fTau(void);

	//establish the size of filter. number of elements contained in a period of the thruster input
	int size_filter (void);

	//do the mean value of parameters in a period. Used as a low pass filter.
	void filter_parameters(int size);

	//normalized error. calculated as the mean absolute velocity error divided by the mean absolute velocity in a period of time (size_filter)
	void mean_norm_error(base::Vector6d _velocity, int size);

	void parameters_estimation(base::VectorXd _thrusterInput, base::Vector6d _velocity);

	void establish_dof(base::VectorXd _thrusterInput, base::Vector6d _velocity);

	base::Vector6d forces_torques (base::VectorXd thrusterInput);

	bool gainsOk;
	bool definedDof;

private:
	//Vector4d m_conventionalParamenters;
	double estimatedVelocity;
	double deltaVelocity;

	base::Vector4d estimatedPhi;
	base::Vector4d estimatedF;
	base::Vector6d gainA;

	//Matrix6x4 gainLambda;
	Eigen::Matrix<double, 6, 4, Eigen::DontAlign> gainLambda;
	base::Vector4d filteredParametersModel;
	base::Vector4d parametersModel;
	double step;
	double fTau;			// frequency of the thruster input in rad/s
	DOFS dof;
	std::queue<base::Vector4d> queueOfParameters;
	std::queue<double> queueOfErrorVelocity;
	std::queue<double> queueOfVelocity;
	double meanErrorVelocity;
	double meanVelocity;
	double normErrorVelocity;

	Eigen::MatrixXd thrusterMatrix;

	int filter_size;
	//int size_filter;

};
}

#endif /* ADAP_PARAMETERS_HPP_ */



















