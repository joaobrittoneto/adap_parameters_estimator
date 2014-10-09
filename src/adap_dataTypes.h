/***************************************************************************/
/*  Data types for an the adaptive parameter estimator of underwater vehicle 	                           */
/*                                                                         */
/* FILE --- uwv_dataTypes.h		                                   */
/*                                                                         */
/* PURPOSE --- Header file for a data types used in the adapt parameters
/ * estimation method for the motion model 					   */
/*                                                                         */
/*  Joao Neto                                             */
/*  joao.neto@dfki.de                                               */
/*  DFKI - BREMEN 2014                                                     */
/***************************************************************************/
#ifndef _ADAP_DATATYPES_H_
#define _ADAP_DATATYPES_H_

#include <base/eigen.h>
#include <base/Eigen.hpp>
#include <math.h>
#include <vector>



namespace adap_parameters_estimator
{
		enum DOFS
			{
				SURGE = 0,
				SWAY  = 1,
				HEAVE = 2,
				ROLL  = 3,
				PITCH = 4,
				YAW   = 5,
				UNINITIALISED = 6,

			};

		struct Direction
			{
				double positive;						// thruster coefficient in positive direction
				double negative;						// thruster coefficient in negative direction
			};

		struct Parameters
			{
				Direction inertiaCoeff[6];			// Coefficients for the diagonal inertial matrix
				Direction linearDampingCoeff[6];
				Direction quadraticDampingCoeff[6];
				base::VectorXd gravityAndBuoyancy;
				base::MatrixXd coriolisCentripetalMatrix;
			};


		//typedef Eigen::Matrix<double, 6, 1, Eigen::DontAlign> Vector6d;
		//typedef Eigen::Matrix<double, 4, 1, Eigen::DontAlign> Vector4d;
		//typedef Eigen::Matrix<double, 4, 4, Eigen::DontAlign> Matrix4d;
		//typedef Eigen::Matrix<double, 6, 4, Eigen::DontAlign> Matrix6x4;

		//<double, 6, 4, Eigen::DontAlign>
}
#endif



















