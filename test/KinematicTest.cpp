/** Standard libraries **/
#include <map>
#include <vector>
#include <string>

/** Library includes **/
#include <threed_odometry/KinematicKDL.hpp>
#include <threed_odometry/KinematicModel.hpp>

#define BOOST_TEST_MODULE KinematicModelTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#ifndef D2R
#define D2R M_PI/180.00 /** Convert degree to radian **/
#endif
#ifndef R2D
#define R2D 180.00/M_PI /** Convert radian to degree **/
#endif

using namespace threed_odometry;

BOOST_AUTO_TEST_CASE( URDFModel)
{

    std::string urdf_file = boost::unit_test::framework::master_test_suite().argv[1]; // Takes path for URDF file from arguments

    /***************************************/
    /** ATTENTION: You need to change the **/
    /** following values to the right ones**/
    /** for your Robot model. Manually in **/
    /** the code.                         **/
    /***************************************/
    std::vector<float> wheel_radius;
    int _RobotTrees, _SlipDoF, _ContactDoF;

    /** VALUES FOR EXOTER **/
    _RobotTrees = 6; _SlipDoF = 3; _ContactDoF = 1;
    wheel_radius.resize(_RobotTrees);
    for(std::vector<float>::iterator it = wheel_radius.begin(); it != wheel_radius.end(); ++it)
    {
        *it = 0.072703; //Wheel radius value
    }
    /*****************************************/
    /** END OF EDITABLE PART FOR PARAMETERS **/
    /*****************************************/

    std::cout<<"********* ROBOT KDL MODEL *********\n";
    KinematicKDL robotKDL(urdf_file, _RobotTrees, _SlipDoF, _ContactDoF, wheel_radius);
    std::cout<<"** URDF FILE: "<<urdf_file<<"\n";
    std::cout<<"** ROBOT MODEL_DOF: "<< robotKDL.MODEL_DOF <<"\n";
    for (std::vector<float>::iterator it = wheel_radius.begin() ; it != wheel_radius.end(); ++it)
        std::cout<<"** ROBOT WHEELS RADIUS: "<< *it <<"\n";

    std::vector<double> joints (robotKDL.MODEL_DOF, 0);
    std::vector<Eigen::Affine3d> fkRobotKDL; // Vector of Forward Kinematics Transformation
    std::vector<base::Matrix6d> fkCovKDL; // Vector of Covariance matrices

    /** Store the values **/
    for (register int i=0; i<static_cast<int>(robotKDL.MODEL_DOF); ++i)
    {
        joints[i] = 0.00;
    }

    std::cout<<"** ROBOT JOINTS VALUES\n";
    for (std::vector<double>::iterator it = joints.begin() ; it != joints.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';

    std::cout<<"** ROBOT MAX CHAIN DoF: "<<robotKDL.getMaxChainDoF()<<"\n";

    robotKDL.fkSolver(joints, fkRobotKDL, fkCovKDL);

    std::cout<<"** fkRobotKDL is of size: "<<fkRobotKDL.size()<<" fkCovKDL of size: "<<fkCovKDL.size()<<"\n";

    for (register unsigned int i=0; i<fkRobotKDL.size(); ++i)
    {
        std::cout<<"fkRobotKDL chain["<<i<<"]\n"<<fkRobotKDL[i].matrix()<<"\n";
    }

    std::cout<<"**\n\n******* ROBOT KDL JACOBIAN *********\n";
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> JKdl;
    JKdl.resize(6*_RobotTrees, robotKDL.MODEL_DOF);

    JKdl = robotKDL.jacobianSolver(joints);

    std::cout<<"** JACOBIAN KDL is of size "<<JKdl.rows()<<" x "<<JKdl.cols()<<"\n"<< JKdl <<"\n\n";
}
