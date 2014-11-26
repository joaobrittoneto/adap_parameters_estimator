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
	~AdapParameters();

	void configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau);

	//Convert the parameters used in the adaptive law for those used in the motion model
	void convetional_parameters (base::VectorXd &estimatedPhi, base::VectorXd &parametersModel);

	// Verify is the gains would lead to stability (gain_a<0, gain_lambda>0)
	void check_gains (void);

	//Construct the diagonal matrix of gain lambda for the specified dof.
	void matrix_lambda (base::Matrix4d &matrixLambda);

	//establish the size of filter. number of elements contained in a period of the thruster input
	int size_filter (void);

	// Simple Moving Average filter.
	void SMA(std::queue<base::VectorXd> &queue, base::VectorXd &filteredValue);

	void parameters_estimation(base::VectorXd &_thrusterInput, base::Vector6d &_velocity, base::Vector4d &estimatedParameters, double &deltaV, double &norm_Error);

	void establish_dof(base::VectorXd _thrusterInput, base::Vector6d _velocity);

	void forces_torques (base::VectorXd &thrusterInput, base::Vector6d &forcesTorques);


	bool gainsOk;
	bool definedDof;


	// Construct a queue of fixed size
	template<typename Type>
	void Queue (int size, Type sample, std::queue<Type> &queue)
	{	// number of parameters use to get the mean value.
		// Fill the queue with n samples
		if (queue.size() < size)
		{	//add a new element in the queue
			queue.push (sample);
		}

		if (queue.size() > size && !queue.empty())
		{	//add a new element in the queue while reduce its size by two till the the queue reach a smaller size
			//remove least element
			queue.pop ();

			//insert new element
			queue.push (sample);

			//remove least element
			queue.pop ();
		}

		if (queue.size() == size && !queue.empty())
		{	//remove least element
			queue.pop();
			//insert new element
			queue.push (sample);
		}
	}


private:

	// configure variables
	base::Vector6d gainA;
	Eigen::Matrix<double, 6, 4, Eigen::DontAlign> gainLambda;
	double step;
	double fTau;			// frequency of the thruster input in rad/s
	DOFS dof;
	Eigen::MatrixXd thrusterMatrix;

	bool reset_values;


};
}

#endif /* ADAP_PARAMETERS_HPP_ */



















