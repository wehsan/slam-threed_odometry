rock_library(threed_odometry
    SOURCES KinematicKDL.cpp
    DEPS_PKGCONFIG base-types base-lib kdl_parser
    HEADERS KinematicKDL.hpp
            MotionModel.hpp
            IIR.hpp)

