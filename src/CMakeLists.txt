rock_library(3d_odometry
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(3d_odometry_bin Main.cpp
    DEPS 3d_odometry)

