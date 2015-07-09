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

	//Configure the library.
	void AdapParameters::configure (Eigen::Matrix<double, 6, 4, Eigen::DontAlign> _gainLambda, Eigen::Matrix<double, 6, 1, Eigen::DontAlign> _gainA, DOFS _dof, double _frequencyTau)
	{
		if (gainA!=_gainA || gainLambda!=_gainLambda || fTau!=_frequencyTau || dof!=_dof)
		{
			gainA = _gainA;
			gainLambda = _gainLambda;
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

	void AdapParameters::update_step(double _sampTime)
	{
		if(step!=_sampTime)
			step = _sampTime;
	}

	// Convert the parameters used in the adaptive method to the conventional parameters use in the model
	void AdapParameters::convetional_parameters(base::VectorXd &estimatedPhi, base::VectorXd &parametersModel)
	{
		//Limit of the inertia (uwv+added mass) in each dof is 10000. ParametersModel(0) can't be less or equal 0 (imply in negative or infinite mass)
		double alpha = estimatedPhi(0);
		if (alpha < 0.00001)
			alpha = 0.00001;

		parametersModel(0) = 1/alpha;
		parametersModel(1) = -estimatedPhi(1)/alpha;
		parametersModel(2) = -estimatedPhi(2)/alpha;
		parametersModel(3) = -estimatedPhi(3)/alpha;

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

	//Make the initial try with the dry-mass or momentum of inertia. Need experimental validation
	base::VectorXd AdapParameters::initial_guess (void)
	{

		double mass	= 100;
		double lx	= 1.1;
		double ly	= 0.8;
		double lz	= 0.4;

		double parameter;
		double momentum;
		base::VectorXd init = base::VectorXd::Zero(4);
		if(mass > 0)
		{
			if(dof < 3)
			{
				parameter = 1/mass;
			}
			// Considering a rectangular parallelepiped (1/12*m*(l1²+l2²))
			// roll
			else if(dof==3)
			{
				momentum = mass*(ly*ly + lz*lz)/12;
				parameter = 1/momentum;
			}
			// pitch
			else if(dof==4)
			{
				momentum = mass*(lx*lx + lz*lz)/12;
				parameter = 1/momentum;
			}
			// yaw
			else if(dof==5)
			{
				momentum = mass*(lx*lx + ly*ly)/12;
				parameter = 1/momentum;
			}
			else
				parameter = 0;
		}
		else
			parameter = 0;

		init[0] = parameter;
		return init;
	}


	// Simple Moving Average filter.
	void AdapParameters::SMA(std::queue<base::VectorXd> &queueInput, base::VectorXd &filteredValue)
	{
		std::queue<base::VectorXd> queue = queueInput;
		filteredValue = base::VectorXd::Zero(filteredValue.size());
		base::VectorXd temp;
		for (int i=0; i < queueInput.size(); i++)
		{
			temp = queue.front();
			queue.pop();
			filteredValue += temp;
		}
		filteredValue /= queueInput.size();
	}


	void AdapParameters::parameters_estimation(base::Vector6d &_forcesTorques, base::Vector6d &_velocity, base::Vector4d &estimatedParameters, double &deltaV, double &norm_Error)
	{

		// states
		static double estimatedVelocity 				= 0;
		static double estimatedAcceleration 			= 0;
		static base::VectorXd estimatedPhi				= base::VectorXd::Zero(4);
		//static base::VectorXd estimatedPhi				= initial_guess(); // Need test
		static base::VectorXd estimatedPhiDot			= base::VectorXd::Zero(4);

		// Aux variables
		double deltaVelocity 							= 0;
		double meanErrorVelocity 						= 0;
		double meanVelocity 							= 0;
		double normErrorVelocity 						= 0;

		base::VectorXd filteredParametersModel			= base::VectorXd::Zero(4);
		base::VectorXd parametersModel					= base::VectorXd::Zero(4);

		base::Vector4d estimatedF						= base::Vector4d::Zero(4);
		base::VectorXd normError						= base::VectorXd::Zero(2);
		//static base::VectorXd deltaPhi					= base::VectorXd::Zero(4);

		// Queue of errors
		//static std::queue<base::VectorXd> queueOfParameters;
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
			//estimatedPhi					= initial_guess();
			estimatedPhiDot					= base::VectorXd::Zero(4);
			estimatedF						= base::Vector4d::Zero(4);
			normError						= base::VectorXd::Zero(2);

//			while (!queueOfParameters.empty())
//			{
//				queueOfParameters.pop();
//			}
			while (!queueOfNormError.empty())
			{
				queueOfNormError.pop();
			}

			reset_values = false;
		}

		static bool first_time = true;
		if(first_time)
			{
				base::Vector4d phi = base::Vector4d::Zero(4);
				phi = {0.01 ,0, 0, 0};
				//estimatedPhi[0]	= 0.01;// phi;
				first_time		= false;
				std::cout << "massIniti " << estimatedPhi[0] << std::endl;
			}


		base::Matrix4d matrixLambda;
		matrix_lambda(matrixLambda);

		// size of the filter, based on the frequency of the force input and in the sample time.
		int size = size_filter();

		///// Calcu states
		// Euler integrator of velocity
		estimatedVelocity = estimatedVelocity + estimatedAcceleration * step;
		// Euler integrator of parameters
		estimatedPhi = estimatedPhi + estimatedPhiDot * step;
		//deltaPhi = estimatedPhiDot * step;

		///// Update states
		// Error of velocity
		deltaVelocity = estimatedVelocity - _velocity(dof);
		// Estimated states
		estimatedF = {_forcesTorques(dof), (estimatedVelocity*fabs(estimatedVelocity)), estimatedVelocity, 1};
		// Adaptive estimator
		estimatedAcceleration = (gainA(dof) * deltaVelocity) + (estimatedPhi.transpose() * estimatedF);
		// Update law
		estimatedPhiDot = -matrixLambda * deltaVelocity * estimatedF;

		// Transform into conventional parameters (inertia, damping, buoyancy)
		convetional_parameters(estimatedPhi, parametersModel);


		//Queue(size, parametersModel, queueOfParameters);
		//SMA(queueOfParameters, filteredParametersModel);
		//estimatedParameters = filteredParametersModel;

		// Mean error
		normError[0]=fabs(deltaVelocity); normError[1]=fabs(double(_velocity(dof)));
		Queue(size, normError, queueOfNormError);
		SMA(queueOfNormError, normError);

		// Output values
		estimatedParameters	=	parametersModel;
		deltaV = deltaVelocity;
		norm_Error = normError[0]/normError[1];

		// Lyapunov function
		//double W = deltaVelocity*deltaVelocity + deltaPhi.transpose() * matrixLambda.inverse() * deltaPhi;
		//norm_Error = W;

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




















