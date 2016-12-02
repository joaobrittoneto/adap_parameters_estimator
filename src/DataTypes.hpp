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
        SURGE = 0,
        SWAY  = 1,
        HEAVE = 2,
        ROLL  = 3,
        PITCH = 4,
        YAW   = 5,
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
