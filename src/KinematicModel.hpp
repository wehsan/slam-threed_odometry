/**@file KinematicModel.hpp
 *
 * @brief Kinematic Model Abstract class
 *
 * This is the Template Abstract class for the Forward kinematics of  Wheel or Leg Robot
 * considering all the DOFs in the Kinematics Chains (Trees). Basically the forward kinematics
 * transform wheel/Leg contact points position with respect to the body_frame considering the
 * joint values (chain configuration). The wheelJacobian or legJacobian works in velocity space
 * and gives information of each wheel/leg contact point to the contribution of the movement
 * (rover angular and linear velocities).
 *
 * It has been designed to use together with the MotionModel Solver (Least-Squares solver)
 * implemented in the MotionModel.hpp file of this library.
 *
 * It should first compute the forward kinematics in order to have afterwards the Robot Jacobian.
 * Neither the forward kinematics nor the Jacobian are provided by this class since this is an
 * abstract class and the kinematics is platform dependent. Your class should make this job.
 *
 *
 * The rational behind is to have a wrapper abstract class to solve the statistical motion model
 * of any complex robot. The goal is for odometry (localization) purpose in a unify way.
 * The real implementation of the Robot Kinematic Model which would inherits from this class can use
 * any method. This means that either analytical or numerical solutions can be internally apply.
 * An example of numerical solution is the KDL library. An example of analytical is your own closed-form.
 * Since forward kinematics is normally quite easy to get an analytical form is recommended.
 * However, if your robot kinematics is complex
 * and for control purpose you have a Tree structure representation for numerical
 * inverse kinematics in KDL/or similar. Perhaps, it is worthy to wrap this kinematic tree
 * implementation with this class.
 *
 * Anyway, the important aspect is to commit the interface. The nomenclature attempts to be generic.
 * An important issue with respect to a conventional kinematics is that a 3D slip vector
 * is modeled as part of the DoF of the kinematics in the form of displacement/rotation
 * vector (x, y and rot theta). Located at the frame between the robot contact point (feet/wheel)
 * and the ground.
 *
 * Further Details at:  P.Muir et. al Kinematic Modeling of Wheeled Mobile Robots
 *                      M. Tarokh et. al Kinematics Modeling and Analyses of Articulated Rovers
 *                      J. Hidalgo et. al Kinematics Modeling of a Hybrid Wheeled-Leg Planetary Rover
 *                      J. Hidalgo, Navigation and Slip Kinematics for High Performance Motion Models
 *
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date May 2013.
 * @version 1.0.
 */



#ifndef ODOMETRY_KINEMATIC_MODEL_HPP
#define ODOMETRY_KINEMATIC_MODEL_HPP

#include <string> // std string library
#include <vector> // std vector library
#include <base/eigen.h> //base eigen of rock

namespace threed_odometry
{
    template <typename _Scalar> class KinematicModel
    {
        protected:
            int RobotTrees, RobotJointDoF, SlipDoF, ContactDoF;

        public:
            /** The maximum number of DoF that the robot has in the Joint Space (Joints + Slip + Contact angle */
            unsigned int MAX_CHAIN_DOF;

            /** The complete number of DoF of the Robot Model. This is compute for the whole Robot Jacobian Matrix */
            unsigned int MODEL_DOF;

        protected:
            /**@brief return @MAX_CHAIN_DOF. Equivalent to the number of columns in the Jacobian.
             */
            inline void setMaxChainDoF(const unsigned int max_chain_dof)
            {
                this->MAX_CHAIN_DOF = max_chain_dof;
            }
            /**@brief return @MODEL_DOF.
             */
            inline void setModelDoF(const unsigned int model_dof)
            {
                this->MODEL_DOF = model_dof;
            }

        public:

            /**@brief Constructor
             **/
            KinematicModel(const int _RobotTrees, const int  _RobotJointDoF, const int _SlipDoF, const int _ContactDoF)
                            :RobotTrees(_RobotTrees), RobotJointDoF(_RobotJointDoF), SlipDoF(_SlipDoF), ContactDoF(_ContactDoF)
            {
                MAX_CHAIN_DOF = -1;
                MODEL_DOF = -1;
            }

            /** @brief return the number of Trees @RobotTrees
             */
            inline int getNumberOfTrees()
            {
                return RobotTrees;
            }
            /** @brief return the complete robot number of Joints @RobotJointDoF
             */
            inline int getRobotJointDoF(void)
            {
                return RobotJointDoF;
            }
            /** @brief return the dimension of the slip model @SlipDoF
             */
            inline int getSlipDoF(void)
            {
                return SlipDoF;
            }
            /**@brief return the dimension of the contact angle model @ContactDoF
             */
            inline int getContactDoF()
            {
                return ContactDoF;
            }
            /**@brief return @MAX_CHAIN_DOF. Equivalent to the number of columns in the Jacobian.
             */
            inline unsigned int getMaxChainDoF()
            {
                return MAX_CHAIN_DOF;
            }
            /**@brief return @MODEL_DOF.
             */
            inline unsigned int getModelDoF()
            {
                return MODEL_DOF;
            }

            /**@brief The name of the Kinematic Model
             */
            virtual std::string name() = 0;

            /**@brief Number of contact points per Tree
             *
             * Candidate points in contact with the groudn. For example in the case of
             * <a href="http://robotik.dfki-bremen.de/en/forschung/robotersysteme/asguard-ii.html">Asguard</a>
             * all the Trees have the same number of contact poins (5).
             * 
             * @param[out] contactPoints the vector of size number of Trees with the number of contact points
             */
            virtual void contactPointsPerTree (std::vector<unsigned int> &contactPoints) = 0;

            /**@brief Assign the current point in contact
             *
             * Point in contact among the number of contact points.
             *
             * @param[in] pointsInContact vector of size number of Trees with the point in contact per each Tree
             */
            virtual void setPointsInContact (const std::vector<int> &pointsInContact) = 0;

            /**@brief Forward Kinematics Solver
             *
             * Compute the Forward Kinematics per each Tree of the robot.
             * Since Forward Kinematics only makes sense in chains and not in Tree. The selected chain per each Tree is that one with the current 
             * point in contact @sa{setPointsInContact}
             *
             * @param[in] positions vector of size @sa{MODEL_DOF} with the joint position values of the complete robot.
             * @param [out] Eigen::Affine3d> &fkRobot, vector of homogeneus transformation of the FK per Tree.
             * @param[out] std::vector<base::Matrix6d> &fkCov, vector of noise covariance associated to each transformation.
             */
            virtual void fkSolver(const std::vector<_Scalar> &positions, std::vector<Eigen::Affine3d> &fkRobot, std::vector<base::Matrix6d> &fkCov) = 0;

            /**@brief Computes the Robot Jacobian
             *
             * Computes the robot Jacobian. Dimension 6*@sa{_RobotTrees} x @sa{MODEL_DOF}
             *
             * @param[in] positions with the current values joints values
             *
             * @return the robot Jacobian matrix.
             */
            virtual Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> jacobianSolver(const std::vector<std::string> &names, const std::vector<_Scalar> &positions) = 0;
    };
}

#endif // ODOMETRY_KINEMATIC_MODEL_HPP

