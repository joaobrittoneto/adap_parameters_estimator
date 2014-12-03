#include <adap_parameters_estimator/adap_parameters.hpp>


//using namespace seabotix_lbv150;
using namespace adap_parameters_estimator;

int main(int argc, char** argv)
{
	base::Vector6d gainA;
	Eigen::Matrix<double, 6, 4, Eigen::DontAlign> gainLambda;
	base::Vector6d tau;
	base::Vector6d realVelocity;
	base::MatrixXd thrusterMatrix;
	//Eigen::Matrix<double, 6, 5, Eigen::DontAlign> thrusterMatrix;

	double step = 0.01;
	double ftau = 1;
	int interactions = 0;

	Parameters param;



	//estimatedPhi <<  0, 0, 0, 0;
	//estimatedF << 0, 0, 0, 0;
	//estimatedVelocity = 0;

	gainA << -0.1, -0.0001, -1, -1, -1, -0.001;

	gainLambda << 	0.0005, 0.5, 0.5, 0.0005,
					0.0001, 1,   1,   0,
					1,      1,   1,   1,
					1,      1,   1,   1,
					1,      1,   1,   1,
					0.02,   0.5, 0.5, 0.0005;



	thrusterMatrix.resize(6,5);

	thrusterMatrix <<	1, 1, 0, 0, 0,
						0, 0, 1, 0, 0,
						0, 0, 0, 1, 1,
						0, 0, 0, 0, 0,
						0, 0, 0, 0, 0,
						-0.21, 0.21, -0.3, 0, 0;



	AdapParameters adapParameters(gainLambda, gainA, thrusterMatrix, SURGE, step, ftau);



	for(double t=0; t<600; t=t+step)
	{

			//realVelocity(SURGE) = 1;
			//tau(SURGE) = 100;
			//TODO descobrir velocidade para frequencia diferente de 1hz
			realVelocity(SURGE) = 0.519*sin( adapParameters.get_fTau()*t-1)+1.53;
			tau(SURGE) = 2*(50*sin(adapParameters.get_fTau()*t)+50);

			realVelocity(SWAY) = 0.2505*sin(adapParameters.get_fTau()*t-1)+0.785;
			tau(SWAY) = (50*sin(adapParameters.get_fTau()*t)+50);

			realVelocity(YAW) = 0.74*sin(adapParameters.get_fTau()*t-0.45)+1.34;
			tau(YAW) = (15*sin(adapParameters.get_fTau()*t)+15);

			adapParameters.parameters_estimation(tau, realVelocity);

			/*
			adapParameters.delta_velocity(realVelocity);
			adapParameters.estimated_state(tau);
			adapParameters.euler_velocity();
			adapParameters.euler_parameters();
			adapParameters.convetional_parameters();
			 */

			adapParameters.filter_parameters(100);

			

			std::cout << std::endl << "===================="<< std::endl;
			//std::cout << std::endl << "Acceleration: " <<adapParameters.estimatedAcceleration  << std::endl << std::endl;
			//std::cout << std::endl << "'plant': " <<adapParameters.estimatedPhi.transpose()*estimatedF  << std::endl << std::endl;
			//std::cout << std::endl << "phi "<<adapParameters.estimatedPhi.transpose() << std::endl << std::endl;
			//std::cout << std::endl <<"EST:"<< adapParameters.get_est_velocity() << std::endl << std::endl;
			//std::cout << std::endl <<"EST²:"<< adapParameters.estimatedVelocity*fabs(adapParameters.estimatedVelocity) << std::endl << std::endl;
			std::cout << std::endl <<"DELTA_v: "<< adapParameters.get_delta_v() << std::endl << std::endl;
			std::cout << std::endl << "Estimated States:  "<< std::endl << adapParameters.get_est_states() << std::endl << std::endl;
			//std::cout << std::endl << adapParameters.tau(SURGE) << std::endl << std::endl;
			std::cout << std::endl << "Filtered Parameters: "<< std::endl <<adapParameters.get_filtered_parameters() << std::endl << std::endl;
			//std::cout << std::endl << "REAL velocity: "<<realVelocity(SURGE) << std::endl << std::endl;
			//std::cout << std::endl <<"Parameters " << std::endl << adapParameters.get_parameters() << std::endl << std::endl;
			//std::cout << std::endl <<"Filtered Parameters " << std::endl << adapParameters.get_filtered_parameters() << std::endl << std::endl;
	}

	//std::cout << std::endl << adapParameters.matrix_lambda(gainLambda, SURGE) << std::endl << std::endl;
	//std::cout << std::endl << adapParameters.get_est_phi() << std::endl << std::endl;
	//adapParameters.convetional_parameters();
	//std::cout << std::endl << adapParameters.get_filtered_parameters() << std::endl << std::endl;

	param.inertiaCoeff[SURGE].positive = adapParameters.get_filtered_parameters()[0];
	param.quadraticDampingCoeff[SURGE].positive = adapParameters.get_filtered_parameters()[1];
	param.linearDampingCoeff[SURGE].positive = adapParameters.get_filtered_parameters()[2];
	param.gravityAndBuoyancy[SURGE] = adapParameters.get_filtered_parameters()[3];

	//std::cout << std::endl << param.inertiaCoeff[SURGE].positive << std::endl << std::endl;
	//std::cout << std::endl << param.quadraticDampingCoeff[SURGE].positive << std::endl << std::endl;
	//std::cout << std::endl << param.linearDampingCoeff[SURGE].positive << std::endl << std::endl;
	//std::cout << std::endl << param.gravityAndBuoyancy[SURGE] << std::endl << std::endl;

	return 0;

}

