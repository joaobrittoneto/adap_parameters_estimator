#ifndef _ADAP_DATATYPES_H_
#define _ADAP_DATATYPES_H_

#include <base/Time.hpp>
#include <base/Eigen.hpp>
#include <math.h>
#include <vector>

namespace adap_parameters_estimator
{
    enum DOF
    {
        LINEAR_X  = 0,
        LINEAR_Y  = 1,
        LINEAR_Z  = 2,
        ANGULAR_X = 3,
        ANGULAR_Y = 4,
        ANGULAR_Z = 5,
    };

    struct OneDOFParameters
    {
        DOF dof;
        base::Time time;
        double inertia;
        double linear_damping;
        double quadratic_damping;
        double buoyancy;
    };
}
#endif
