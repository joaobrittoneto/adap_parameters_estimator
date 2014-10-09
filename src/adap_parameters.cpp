/***************************************************************************/
/*  Adaptive Parameters Identification Library for AUV                     */
/*                                                                         */
/* FILE --- adap_parameters.cpp	                                           */
/*                                                                         */
/* PURPOSE --- Source file for a Adaptive Parameters Identification to be  */
/* used in a motion model of a underwater vehicle.	                       */
/*             Based on D. A. Smallwood(2003)                              */
/*                                                                         */
/*  João Da Costa Britto Neto                                             */
/*  joao.neto@dfki.de                                               */
/*  DFKI - BREMEN 2014                                                     */
/***************************************************************************/

#include "adap_parameters.hpp"


namespace adap_parameters_estimator{

//TODO IMPLEMENTAR UM SO CONSTRUTOR E POR E USAR O DEFAULT DO OROGEN PARA SETAR O DOF


	AdapParameters::AdapParameters(Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau)
		{

			gainA = _gainA;
			gainLambda = _gainLambda;
			step = _sampTime;
			fTau = _frequencyTau;
			thrusterMatrix = _thrusterMatrix;
			dof = _dof;

			if (dof == UNINITIALISED)
				definedDof = false;
			else
				definedDof = true;
			gainsOk = true;
			check_gains();

			filter_size = 0;

			estimatedVelocity = 0;
			deltaVelocity = 0;
			normErrorVelocity = 0;
			meanErrorVelocity = 0;
			meanVelocity = 0;
			estimatedPhi << 0, 0, 0, 0;
			estimatedF << 0, 0, 0, 0;
			parametersModel << 0, 0, 0, 0;
			filteredParametersModel << 0, 0, 0, 0;
		}

/*	AdapParameters::AdapParameters()
	{
		estimatedVelocity = 0;
		deltaVelocity = 0;
		estimatedPhi << 0, 0, 0, 0;
		estimatedF << 0, 0, 0, 0;
		parametersModel << 0, 0, 0, 0;
		filteredParametersModel << 0, 0, 0, 0;
	}
*/

	AdapParameters::~AdapParameters()
	{
	}

	//Configure the library.
	void AdapParameters::configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau)
	{
		gainA = _gainA;
		gainLambda = _gainLambda;
		step = _sampTime;
		fTau = _frequencyTau;
		thrusterMatrix = _thrusterMatrix;
		dof = _dof;

		if (dof == UNINITIALISED)
			definedDof = false;
		else
			definedDof = true;
		gainsOk = true;
		check_gains();

		filter_size = 0;
		estimatedVelocity = 0;
		deltaVelocity = 0;
		normErrorVelocity = 0;
		meanErrorVelocity = 0;
		meanVelocity = 0;
		estimatedPhi << 0, 0, 0, 0;
		estimatedF << 0, 0, 0, 0;
		parametersModel << 0, 0, 0, 0;
		filteredParametersModel << 0, 0, 0, 0;

		while (!queueOfParameters.empty())
		{
			queueOfParameters.pop();
		}
		while (!queueOfErrorVelocity.empty())
		{
			queueOfErrorVelocity.pop();
		}
		while (!queueOfVelocity.empty())
		{
			queueOfVelocity.pop();
		}

	}

	//Construct the vector of non-linear states
	void AdapParameters::estimated_state (base::Vector6d _tau)
	{
		estimatedF = {_tau(dof), (estimatedVelocity*fabs(estimatedVelocity)), estimatedVelocity, 1};

	}

	void AdapParameters::delta_velocity(base::Vector6d _velocity)
	{
		deltaVelocity = estimatedVelocity - _velocity(dof);
		//normDeltaVelocity = 100*fabs(deltaVelocity)/fabs(_velocity(dof));
	}

	// Adaptive estimator
	double AdapParameters::adaptive_estimator (void)
	{
		double estimatedAcceleration;

		estimatedAcceleration = (gainA(dof) * deltaVelocity) + (estimatedPhi.transpose() * estimatedF);

		return estimatedAcceleration;
	}

	// Update law of the parameters
	base::Vector4d AdapParameters::update_law(void)
	{
		base::Vector4d estimatedPhiDot;
		base::Matrix4d matrixLambda = matrix_lambda();

		estimatedPhiDot = -matrixLambda * deltaVelocity * estimatedF;

		return estimatedPhiDot;
	}

	// Convert the parameters used in the adaptive method to the conventional parameters use in the model
	void AdapParameters::convetional_parameters(void)
	{
		//Limit of the inertia (uwv+added mass) in each dof is 10000. ParametersModel(0) can't be less or equal 0 (imply in negative or infinite mass)
		if (estimatedPhi(0) >= 0.00001)
			{
			parametersModel(0) = 1/estimatedPhi(0);
			parametersModel(1) = -estimatedPhi(1)/estimatedPhi(0);
			parametersModel(2) = -estimatedPhi(2)/estimatedPhi(0);
			parametersModel(3) = -estimatedPhi(3)/estimatedPhi(0);
			}
		else
			{
			parametersModel(0) = 1;
			parametersModel(1) = 0;
			parametersModel(2) = 0;
			parametersModel(3) = 0;
			}

	}

	// Verify if the gains have the right signal
	void AdapParameters::check_gains(void)
	{	for (int i=0; i<6; i++)
				{
					if(gainA(i) > 0)
						{
						std::cout << std::endl << "The gain_a must be a negative value" << std::endl;
						gainsOk = false;
						}
					else
					{
						for (int j=0; j<4; j++)
						{
							if (gainLambda(i,j) < 0)
								{
								std::cout << std::endl << "The gain_lambda must be a positive value" << std::endl;
								gainsOk = false;
								}


						}
					}
				}
		}

	//Euler method for establish the velocity in the online method
	void AdapParameters::euler_velocity (void)
	{
		double estimatedAccelaration;

		estimatedAccelaration = adaptive_estimator();

		estimatedVelocity = estimatedVelocity + estimatedAccelaration * step;

	}

	//Euler method for establish the parameters in the online method
	void AdapParameters::euler_parameters (void)
	{
		base::Vector4d parametersDot;

		parametersDot = update_law();

		estimatedPhi = estimatedPhi + parametersDot * step;

	}

	//get the diagonal matrix of gain lambda
	base::Matrix4d AdapParameters::matrix_lambda(void)
	{
		base::Vector4d vectorLambda = gainLambda.row(dof);
		base::Matrix4d matrixLambda = vectorLambda.asDiagonal();
		//base::Matrix4d matrixLambda = base::Matrix4d::Zero();
		//matrixLambda[0][0] = gainLambda[dof][0];
		//matrixLambda[1][1] = gainLambda[dof][1];
		//matrixLambda[2][2] = gainLambda[dof][2];
		//matrixLambda[3][3] = gainLambda[dof][3];
		return matrixLambda;
	}

	//get the parameters used in the model
	base::Vector4d AdapParameters::get_parameters (void)
	{
		return parametersModel;
	}

	//return the filtered conventional parameters;
	base::Vector4d AdapParameters::get_filtered_parameters (void)
	{
		return filteredParametersModel;
	}

	//get the error of velocity
	double AdapParameters::get_delta_v (void)
	{
		return deltaVelocity;
	}

	//get the error of velocity
	double AdapParameters::get_norm_delta_v (void)
	{
		return normErrorVelocity;
	}

	//get the nonliner vector of states
	base::Vector4d AdapParameters::get_est_states(void)
	{
		return estimatedF;
	}

	//get the estimated parameters
	base::Vector4d AdapParameters::get_est_phi(void)
	{
		return estimatedPhi;
	}

	//get the estimated velocity
	double AdapParameters::get_est_velocity(void)
	{
		return estimatedVelocity;
	}

	//return frequency of thruster input signal
	double AdapParameters::get_fTau(void)
	{
		return fTau;
	}

	//
	void AdapParameters::filter_parameters(int size)
	{
		// number of parameters use to get the mean value.
		// Fill the queue with n parametersModel
		if (filter_size < size)
		{	//add a new element in the queue
			queueOfParameters.push (parametersModel);
			//compute the new mean value with the new element
			filteredParametersModel = (( filteredParametersModel * (queueOfParameters.size()-1))
									+ queueOfParameters.back()) / queueOfParameters.size();
		}

		if (filter_size > size && !queueOfParameters.empty())
		{	//add a new element in the queue while reduce its size by two till the the queue reach a smaller size

			//remove the influence of the least element of the queue
			filteredParametersModel = (( filteredParametersModel * queueOfParameters.size())
									- queueOfParameters.front()) / (queueOfParameters.size()-1);
			//remove least element
			queueOfParameters.pop ();

			//insert new element
			queueOfParameters.push (parametersModel);
			//compute new mean value with the influence of new element
			filteredParametersModel =  ( (filteredParametersModel * (queueOfParameters.size()-1))
									+ queueOfParameters.back() ) / (queueOfParameters.size());

			//remove the influence of the least element of the queue
			filteredParametersModel = (( filteredParametersModel * queueOfParameters.size())
									- queueOfParameters.front()) / (queueOfParameters.size()-1);
			//remove least element
			queueOfParameters.pop ();

		}

		if (filter_size == size && !queueOfParameters.empty())
		{
			//remove the influence of the least element of the queue
			filteredParametersModel =  ( (filteredParametersModel * queueOfParameters.size())
									- queueOfParameters.front() ) / (queueOfParameters.size()-1);

			//remove least element
			queueOfParameters.pop();
			//insert new element
			queueOfParameters.push (parametersModel);

			//compute new mean value with the influence of new element
			filteredParametersModel =  ( (filteredParametersModel * (queueOfParameters.size()-1))
									+ queueOfParameters.back() ) / (queueOfParameters.size());
		}
	}


	int AdapParameters::size_filter(void)
	{	// number of parameters = T/step; T: period of input signal (thruster forces); step: sample time
		// T = 2pi/f with f in rad/s.
		int n;
		if (fTau > 0)
			{n = (2*M_PI/fTau) / step;}
		else
			{n = 30;} //number of elements in the filter

		if (filter_size < n)
		{
			filter_size ++;
		}

		if (filter_size > n)
		{
			filter_size --;
		}

		return n;
	}


	void AdapParameters::mean_norm_error(base::Vector6d _velocity, int size)
	{	// mean(fabs(deltaV))/mean(fabs(velocity)) in a periody (T) with n samples

		// Fill the queue with n samples
		if (filter_size < size)
		{
			//Error velocity
			//add a new element in the queue
			queueOfErrorVelocity.push (fabs(deltaVelocity));
			//compute the new mean value with the new element
			meanErrorVelocity = (( meanErrorVelocity * (queueOfErrorVelocity.size()-1))
									+ queueOfErrorVelocity.back()) / queueOfErrorVelocity.size();
			//////////////////////////////////////////
			//Velocity
			//add a new element in the queue
			double vdof = _velocity(dof);
			queueOfVelocity.push (fabs(vdof));
			//compute the new mean value with the new element
			meanVelocity = (( meanVelocity * (queueOfVelocity.size()-1))
									+ queueOfVelocity.back()) / queueOfVelocity.size();


		}

		if (filter_size > size && !queueOfErrorVelocity.empty())
		{	//add a new element in the queue while reduce its size by two till the the queue reach a smaller size

			//delta velocity
			//remove the influence of the least element of the queue
			meanErrorVelocity = (( meanErrorVelocity * queueOfErrorVelocity.size())
									- queueOfErrorVelocity.front()) / (queueOfErrorVelocity.size()-1);
			//remove least element
			queueOfErrorVelocity.pop ();

			//insert new element
			queueOfErrorVelocity.push (fabs(deltaVelocity));
			//compute new mean value with the influence of new element
			meanErrorVelocity =  ( (meanErrorVelocity * (queueOfErrorVelocity.size()-1))
									+ queueOfErrorVelocity.back() ) / (queueOfErrorVelocity.size());

			//remove the influence of the least element of the queue
			meanErrorVelocity = (( meanErrorVelocity * queueOfErrorVelocity.size())
									- queueOfErrorVelocity.front()) / (queueOfErrorVelocity.size()-1);
			//remove least element
			queueOfErrorVelocity.pop ();

			/////////////////////////////////////////
			//Velocity
			//remove the influence of the least element of the queue
			meanVelocity = (( meanVelocity * queueOfVelocity.size())
									- queueOfVelocity.front()) / (queueOfVelocity.size()-1);
			//remove least element
			queueOfVelocity.pop ();

			//insert new element
			double vdof = _velocity(dof);
			queueOfVelocity.push (fabs(vdof));
			//compute new mean value with the influence of new element
			meanVelocity =  ( (meanVelocity * (queueOfVelocity.size()-1))
									+ queueOfVelocity.back() ) / (queueOfVelocity.size());

			//remove the influence of the least element of the queue
			meanVelocity = (( meanVelocity * queueOfVelocity.size())
									- queueOfVelocity.front()) / (queueOfVelocity.size()-1);
			//remove least element
			queueOfVelocity.pop ();
		}

		if (filter_size == size && !queueOfErrorVelocity.empty())
		{	//error velocity
			//remove the influence of the least element of the queue
			meanErrorVelocity =  ( (meanErrorVelocity * queueOfErrorVelocity.size())
									- queueOfErrorVelocity.front() ) / (queueOfErrorVelocity.size()-1);

			//remove least element
			queueOfErrorVelocity.pop();
			//insert new element
			queueOfErrorVelocity.push (fabs(deltaVelocity));

			//compute new mean value with the influence of new element
			meanErrorVelocity =  ( (meanErrorVelocity * (queueOfErrorVelocity.size()-1))
									+ queueOfErrorVelocity.back() ) / (queueOfErrorVelocity.size());

			////////////////////////////////////////////////
			//Velocity
			//remove the influence of the least element of the queue
			meanVelocity =  ( (meanVelocity * queueOfVelocity.size())
									- queueOfVelocity.front() ) / (queueOfVelocity.size()-1);

			//remove least element
			queueOfVelocity.pop();
			//insert new element
			double vdof = _velocity(dof);
			queueOfVelocity.push (fabs(vdof));

			//compute new mean value with the influence of new element
			meanVelocity =  ( (meanVelocity * (queueOfVelocity.size()-1))
									+ queueOfVelocity.back() ) / (queueOfVelocity.size());

		 }

		if (meanVelocity != 0)
		{
			normErrorVelocity = meanErrorVelocity/meanVelocity;
		}

	}


	void AdapParameters::parameters_estimation(base::VectorXd _thrusterInput, base::Vector6d _velocity)
	{
		base::Vector6d forcesTorques = forces_torques (_thrusterInput);
		delta_velocity(_velocity);
		estimated_state(forcesTorques);
		euler_velocity();
		euler_parameters();
		convetional_parameters();
		int n = size_filter();
		filter_parameters(n);
		mean_norm_error(_velocity, n);
	}


	void AdapParameters::establish_dof(base::VectorXd _thrusterInput, base::Vector6d _velocity)
	{
		base::Vector6d Tau = forces_torques (_thrusterInput);
		int activeDof = 0; // number of dofs active
		DOFS Dof[7] = {SURGE, SWAY, HEAVE, ROLL, PITCH, YAW, UNINITIALISED};
		int establishDof = 6;

		for (int i=0; i<6; i++)
		{
			double velo = _velocity(i);
			double tau = Tau(i);

			if (fabs(velo) >= 0.001 && fabs(tau) >= 0.001)
			{	activeDof++;
				establishDof = i;
			}

		}

		if (activeDof > 1)
		{
			dof = UNINITIALISED;
			std::cout << std::endl << "Two or more Degree of Freedom working at same time: " << std::endl;
		}

		if (activeDof == 0)
		{
			dof = UNINITIALISED;
			std::cout << std::endl << "Degree of Freedom not identified: " << std::endl;
		}

		else dof = Dof[establishDof];
	}


	base::Vector6d AdapParameters::forces_torques (base::VectorXd thrusterInput)
	{	base::Vector6d forcesTorques = base::Vector6d::Zero();
		//if (thrusterMatrix.size()/6 == thrusterInput.size()) ////In case the input is the forces applied for each thruster
			//forcesTorques = thrusterMatrix * thrusterInput;  //In case the input is the forces applied for each thruster
		forcesTorques = thrusterInput; // In case the input is the forces and torques applied direct to the auv
		//else //In case the input is the forces applied for each thruster
		//	std::cout << std::endl << "Thruster Matrix need be compatible with number of thrusters. " << std::endl;
		return forcesTorques;
	}
}




















