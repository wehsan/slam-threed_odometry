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


#define EIGEN_NO_AUTOMATIC_RESIZING //To avoid automatic resizing of Dynamic Matrices
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO //Zero as default value for Jacobian matrices

/* Base includes */
#include <assert.h> /** Assert **/
#include <sstream>
#include <string> /** This should be already in sstream **/

/** Rock Base Library **/
#include <base/Time.hpp>
#include <base/Eigen.hpp>
#include <base/Logging.hpp>

#include <Eigen/Core>

/* Odometry library */
#include <threed_odometry/KinematicModel.hpp> /** For the Odometry class to inherit **/

/* KDL Parser */
#include <kdl_parser/kdl_parser.hpp>

/** KDL library **/
#include <kdl/frames_io.hpp> /** Defines KDL frames **/
#include <kdl/chainfksolver.hpp> /** Solve forward kinematics **/
#include <kdl/chainfksolverpos_recursive.hpp> /** Position of the end-effector **/
#include <kdl/chainjnttojacsolver.hpp> /** For Jacobian matrix **/


/* URDF */
#include "urdf_parser/urdf_parser.h"
#include <urdf_model/model.h>

namespace threed_odometry
{

    class KinematicKDL : public ::threed_odometry::KinematicModel<double>
    {

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW //Structures having Eigen members

    protected:
        KDL::Tree tree; /** There are as many chains as number of _RobotTrees **/
        std::string model_name; /** Name of the model class **/
        std::vector<float> wheel_radius; /** Number and wheel radius **/

        /** Wheel Jacobian Matrix  (passive joint for Rear wheels/column of zeros for Front wheels, wheel rotation, (3 x 1)slip vector and contact_angle **/
        std::vector< Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> , Eigen::aligned_allocator < Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> > > vector_jacobians;
        std::vector<int>  contact_idx; /** Foot index making the motion (0,...,n) of the wheelidx wheel **/

    public:

        KinematicKDL(const int _RobotTrees, const int _RobotJointDoF,
                        const int _SlipDoF, const int _ContactDoF,
                        std::string &urdf_file, std::vector<float> &wheel_radius);
        ~KinematicKDL();

        std::string name();
        void contactPointsPerTree (std::vector<unsigned int> &contactPoints);
        void setPointsInContact (const std::vector<int> &pointsInContact);
        void fkBody2ContactPointt(const int chainIdx, const std::vector<double> &positions, Eigen::Affine3d &fkTrans, base::Matrix6d &fkCov);
        void fkSolver(const std::vector<double> &positions, std::vector<Eigen::Affine3d> &fkRobot, std::vector<base::Matrix6d> &fkCov);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> jacobianSolver(const std::vector<double> &positions);

    };

    static void printLink(const KDL::SegmentMap::const_iterator& link, const std::string& prefix)
    {
        //std::cout << prefix << "- Segment " << link->second.segment.getName() << " has " << link->second.children.size() << " children" << std::endl;
        for (unsigned int i=0; i<link->second.children.size(); ++i)
            printLink(link->second.children[i], prefix + "  ");
    };

    static std::string getContactPoint (const KDL::SegmentMap::const_iterator& link, const int childrenIdx)
    {
        if (link->second.children.size() == 0)
            return link->second.segment.getName();

        //std::cout<<link->second.segment.getName() << " has "<<link->second.children.size()<<" children\n";
        return getContactPoint(link->second.children[childrenIdx], 0);
    };

};

#endif
