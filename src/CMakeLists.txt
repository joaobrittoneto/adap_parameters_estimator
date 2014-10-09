rock_library(adap_parameters_estimator
    SOURCES adap_parameters.cpp
    HEADERS adap_parameters.hpp adap_dataTypes.h
    DEPS_PKGCONFIG base-types)
#rock_executable(adap_parameters_estimator_bin Main.cpp
