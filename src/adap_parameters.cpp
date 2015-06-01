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
		check_gains();
		reset_values = true;
	}

	AdapParameters::AdapParameters(Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, DOFS _dof, double _sampTime, double _frequencyTau)
	{
		gainA = _gainA;
		gainLambda = _gainLambda;
		step = _sampTime;
		fTau = _frequencyTau;
		dof = _dof;

		if (dof == UNINITIALISED)
			definedDof = false;
		else
			definedDof = true;
		check_gains();
		reset_values = true;
	}



	AdapParameters::~AdapParameters()
	{
	}

	//Configure the library.
	void AdapParameters::configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, Eigen::MatrixXd _thrusterMatrix, DOFS _dof, double _sampTime, double _frequencyTau)
	{
		if (gainA!=_gainA || gainLambda!=_gainLambda || step!=_sampTime || fTau!=_frequencyTau || thrusterMatrix!=_thrusterMatrix || dof!=_dof)
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
			check_gains();
			reset_values = true;
		}
	}

	//Configure the library.
	void AdapParameters::configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, DOFS _dof, double _sampTime, double _frequencyTau)
	{
		if (gainA!=_gainA || gainLambda!=_gainLambda || step!=_sampTime || fTau!=_frequencyTau || dof!=_dof)
		{
			gainA = _gainA;
			gainLambda = _gainLambda;
			step = _sampTime;
			fTau = _frequencyTau;
			dof = _dof;

			if (dof == UNINITIALISED)
				definedDof = false;
			else
				definedDof = true;
			check_gains();
			reset_values = true;
		}
	}

	// Convert the parameters used in the adaptive method to the conventional parameters use in the model
	void AdapParameters::convetional_parameters(base::VectorXd &estimatedPhi, base::VectorXd &parametersModel)
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
	{	gainsOk = true;
		for (int i=0; i<6; i++)
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

	//get the diagonal matrix of gain lambda
	void AdapParameters::matrix_lambda(base::Matrix4d &matrixLambda)
	{
		base::Vector4d vectorLambda = gainLambda.row(dof);
		matrixLambda = vectorLambda.asDiagonal();
	}

	// Size of queue, for the SMA filter, based on the number of samples in a period
	int AdapParameters::size_filter(void)
	{	// number of parameters = T/step; T: period of input signal (thruster forces); step: sample time
		// T = 2pi/f with f in rad/s.
		int n;
		if (fTau > 0)
			{n = (2*M_PI/fTau) / step;}
		else
			{n = 30;} //number of elements in the filter

		return n;
	}

	// Simple Moving Average filter.
	void AdapParameters::SMA(std::queue<base::VectorXd> &queue, base::VectorXd &filteredValue)
	{
		filteredValue = base::VectorXd::Zero(filteredValue.size());
		base::VectorXd temp;
		for (int i=0; i < queue.size(); i++)
		{
			temp = queue.front();
			queue.pop();
			filteredValue += temp;
			queue.push(temp);
		}
		filteredValue /= queue.size();
	}


	void AdapParameters::parameters_estimation(base::Vector6d &_forcesTorques, base::Vector6d &_velocity, base::Vector4d &estimatedParameters, double &deltaV, double &norm_Error)
	{

		static double deltaVelocity 					= 0;
		static double estimatedVelocity 				= 0;
		static double estimatedAcceleration 			= 0;
		static double meanErrorVelocity 				= 0;
		static double meanVelocity 						= 0;
		static double normErrorVelocity 				= 0;

		static base::VectorXd filteredParametersModel	= base::VectorXd::Zero(4);
		static base::VectorXd parametersModel			= base::VectorXd::Zero(4);
		static base::VectorXd estimatedPhi				= base::VectorXd::Zero(4);
		static base::VectorXd estimatedPhiDot			= base::VectorXd::Zero(4);
		static base::Vector4d estimatedF				= base::Vector4d::Zero(4);
		static base::VectorXd normError					= base::VectorXd::Zero(2);

		static std::queue<base::VectorXd> queueOfParameters;
		static std::queue<base::VectorXd> queueOfNormError;

		if (reset_values == true)
		{
			deltaVelocity 					= 0;
			estimatedVelocity 				= 0;
			estimatedAcceleration 			= 0;
			meanErrorVelocity 				= 0;
			meanVelocity 					= 0;
			normErrorVelocity 				= 0;

			filteredParametersModel			= base::VectorXd::Zero(4);
			parametersModel					= base::VectorXd::Zero(4);
			estimatedPhi					= base::VectorXd::Zero(4);
			estimatedPhiDot					= base::VectorXd::Zero(4);
			estimatedF						= base::Vector4d::Zero(4);
			normError						= base::VectorXd::Zero(2);

			while (!queueOfParameters.empty())
			{
				queueOfParameters.pop();
			}
			while (!queueOfNormError.empty())
			{
				queueOfNormError.pop();
			}

			reset_values = false;
		}


		base::Matrix4d matrixLambda;
		matrix_lambda(matrixLambda);

		// size of the filter, based on the frequency of the force input and in the sample time.
		int size = size_filter();

		// Error of velocity
		deltaVelocity = estimatedVelocity - _velocity(dof);
		// Estimated states
		estimatedF = {_forcesTorques(dof), (estimatedVelocity*fabs(estimatedVelocity)), estimatedVelocity, 1};
		// Adaptive estimator
		estimatedAcceleration = (gainA(dof) * deltaVelocity) + (estimatedPhi.transpose() * estimatedF);
		// Euler integrator of velocity
		estimatedVelocity = estimatedVelocity + estimatedAcceleration * step;
		// Update law
		estimatedPhiDot = -matrixLambda * deltaVelocity * estimatedF;
		// Euler integrator of parameters
		estimatedPhi = estimatedPhi + estimatedPhiDot * step;
		// Transform into conventional parameters (inertia, damping, buoyancy)
		convetional_parameters(estimatedPhi, parametersModel);


		Queue(size, parametersModel, queueOfParameters);
		normError[0]=fabs(deltaVelocity); normError[1]=fabs(double(_velocity(dof)));
		Queue(size, normError, queueOfNormError);

		SMA(queueOfParameters, filteredParametersModel);
		SMA(queueOfNormError, normError);


		estimatedParameters = filteredParametersModel;
		//estimatedParameters	=	parametersModel;
		deltaV = deltaVelocity;
		norm_Error = normError[0]/normError[1];

	}

	// Verify method (Mean values of forces and velocity?)
	void AdapParameters::establish_dof(base::Vector6d _forcesTorques, base::Vector6d _velocity)
	{

		int activeDof = 0; // number of dofs active
		DOFS Dof[7] = {SURGE, SWAY, HEAVE, ROLL, PITCH, YAW, UNINITIALISED};
		int establishDof = 6;

		for (int i=0; i<6; i++)
		{
			double velo = _velocity(i);
			double tau = _forcesTorques(i);

			if (fabs(velo) >= 0.1 && fabs(tau) >= 0.1)
			{	activeDof++;
				establishDof = i;
			}

		}

		if (activeDof > 1)
		{
			dof = UNINITIALISED;
			std::cout << std::endl << "Two or more Degree of Freedom working at same time: " << std::endl;
		}

		else if (activeDof == 0)
		{
			dof = UNINITIALISED;
			std::cout << std::endl << "Degree of Freedom not identified: " << std::endl;
		}

		else dof = Dof[establishDof];
		dof = Dof[0]; //TODO remove this line after verify method
	}


	void AdapParameters::forces_torques (base::VectorXd &thrusterInput, base::Vector6d &forcesTorques)
	{	//if (thrusterMatrix.size()/6 == thrusterInput.size()) ////In case the input are the forces applied for each thruster
			forcesTorques = thrusterMatrix * thrusterInput;  //In case the input are the forces applied for each thruster
		//else //In case the input are the forces applied for each thruster
		//	std::cout << std::endl << "Thruster Matrix need be compatible with number of thrusters. " << std::endl;
		//forcesTorques = thrusterInput; // In case the input are the forces and torques applied direct to the auv
	}
}




















