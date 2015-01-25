/**\file KinematicKDL.hpp
*
* This class has the primitive methods for the Forward kinematics and Wheel Jacobian.
* First computes the forward kinematics in order to have afterwards the wheel Jacobians using KDL.
*
* It read a URDF file to create KDL trees using the kdl_parse.
*
* @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
* @date October 2013.
* @version 1.0.
*/


#ifndef ODOMETRY_KINEMATIC_KDL_HPP
#define ODOMETRY_KINEMATIC_KDL_HPP

/* Base includes */
#include <assert.h> /** Assert **/
#include <sstream>
#include <string> /** This should be already in sstream **/

/** Rock Base Library **/
#include <base/Time.hpp>
#include <base/Eigen.hpp>
#include <base/Logging.hpp>

#include <Eigen/Core>

/* KDL Parser */
#include <kdl_parser/kdl_parser.hpp>

/** KDL library **/
#include <kdl/frames_io.hpp> /** Defines KDL frames **/
#include <kdl/treefksolver.hpp> /** Solve forward kinematics **/
#include <kdl/treefksolverpos_recursive.hpp> /** Solve forward kinematics **/
#include <kdl/chainjnttojacsolver.hpp> /** For Jacobian matrix **/
#include <kdl/treejnttojacsolver.hpp> /** For Jacobian matrix **/


/* URDF */
#include <urdf_parser/urdf_parser.h>
#include <urdf_model/model.h>

namespace threed_odometry
{

    class KinematicKDL
    {

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW //Structures having Eigen members

    protected:
        KDL::Tree tree; /** There are as many chains as number of _RobotTrees **/
        std::string model_name; /** Name of the model class **/
        std::vector<std::string> contact_points; /** Number and segment names of contact points **/
        int number_robot_joints, number_slip_joints, number_contact_joints;

    public:
        unsigned int model_dof;


    public:
        KinematicKDL (const std::string &urdf_file, const std::vector<std::string> &contact_points, 
                const int _number_robot_joints, const int _number_slip_joints, const int _number_contact_joints);
        ~KinematicKDL();

        std::string getName();
        void fkSolver(const std::vector<double> &positions, std::vector<Eigen::Affine3d> &fkTrans, std::vector<base::Matrix6d> &fkCov);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobianSolver(const std::vector<std::string> &names, const std::vector<double> &positions);
        void unpackJoints(const std::vector<std::string>& joint_names, const std::vector<double>& joint_positions,
                                const std::vector<std::string>& involved_joints, KDL::JntArray& joint_array);
        void organizeJacobian(const int chainidx, const std::vector<std::string> &joint_names, const std::vector<std::string> &involved_joints,
                                    const Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> &jacobian, Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> &J);

    };

    static void printLink(const KDL::SegmentMap::const_iterator& link, const std::string& prefix)
    {
        std::cout << prefix << "- Segment " << link->second.segment.getName() << " [Joint "<< link->second.segment.getJoint().getName()<<
            "] has " << link->second.children.size() << " children" << std::endl;
        for (unsigned int i=0; i<link->second.children.size(); ++i)
            printLink(link->second.children[i], prefix + "  ");
    };

    static std::string getContactPoint (const KDL::SegmentMap::const_iterator& link, const int childrenIdx)
    {
        if (link->second.children.size() == 0)
            return link->second.segment.getName();

        std::cout<<link->second.segment.getName() << " has "<<link->second.children.size()<<" children\n";
        return getContactPoint(link->second.children[childrenIdx], 0);
    };

};

#endif
