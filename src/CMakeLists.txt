rock_library(adap_parameters_estimator
    SOURCES AdapParameters.cpp
    HEADERS AdapParameters.hpp DataTypes.hpp
    DEPS_PKGCONFIG base-types)
#rock_executable(adap_parameters_estimator_bin Main.cpp
