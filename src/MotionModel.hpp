/**\file MotionModel.hpp
 *
 * This class provided a Motion Model solver for any mobile robot (wheel or leg)
 * with articulated joints.
 *
 * The solver is based on Weighted Least-Squares as minimizing the error to estimate the
 * resulting motion from a robot Jacobian. Robot Jacobian is understood as the sparse Jacobian
 * matrix containing one Jacobian per each Tree of the Robot.
 * Therefore, it is a minimizing problem of a non-linear system which has been previously linearized
 * around the working point (robot current joint position values). This is represented in the
 * current Jacobian matrix of the robot.
 * This class uses the KinematicModel abstract class (in this library) to access the robot Jacobian
 * and compute the statistical motion model (estimated velocities/navigation quantities).
 *
 * Further Details at:  P.Muir et. al Kinematic Modeling of Wheeled Mobile Robots
 *                      M. Tarokh et. al Kinematics Modeling and Analyses of Articulated Rovers
 *                      J. Hidalgo et. al Kinematics Modeling of a Hybrid Wheeled-Leg Planetary Rover
 *                      J. Hidalgo, Navigation and Slip Kinematics for High Performance Motion Models
 *
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date June 2013.
 * @version 1.0.
 */


#ifndef ODOMETRY_MOTION_MODEL_HPP
#define ODOMETRY_MOTION_MODEL_HPP

#include <iostream> /** for std::cout TO-DO REMOVE **/
#include <vector> /** For std vector **/
#include <boost/shared_ptr.hpp> /** For std shared pointer **/
#include <base/logging.h> /** Log message **/
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */
#include <Eigen/Core> /** Core methods of Eigen implementation **/
#include <Eigen/Dense> /** for the algebra and transformation matrices and accessing Matrixblock and corner among others**/
#include <Eigen/Cholesky> /** For the Cholesky decomposition **/

#define DEBUG_PRINTS_ODOMETRY_MOTION_MODEL 1 //TO-DO: Remove this. Only for testing (master branch) purpose

namespace threed_odometry
{

    /**@class MotionModel
     *
     * Motion Model solver class
     *
     * @param _Scalar: is the typename of your desired implementation (float, double, etc..)
     * @param _RobotTrees: is the number of independent Trees connected to your desired Body Center.
     *              Trees is understood as connection of kinematics chains. As an example:
     *              <a href="http://robotik.dfki-bremen.de/en/forschung/robotersysteme/asguard-ii.html">Asguard</a>
     *              hybrid wheels model as Trees. Therefore Asguard has 4 Trees. Each Tree has 5 open
     *              kinematics chains, one per each foot which are potential points in contact with
     *              the ground.
     * @param _RobotJointDoF: Complete number of DoF in the Joint Space of your robot.This is every joint
     *              passive or active that your robot (or the part of your robot related to odometry)
     *              has. This is required for the Jacobian in a general form. Therefore, this would
     *              be part of the number of columns in the Kinematic Model Jacobian.
     * @param _SlipDoF: The number of DoF than you want to model the slip of the contact point. If you
     *           are not interested to model slip velocity/displacement set it to zero. Otherwise,
     *           a common slip model has 3DoF (in X an Y direction and Z rotation)
     * @param _ContactDoF: This is the angle of contact between the ground and the contact point. If this
     *              is not interesting for you model (i.e. the robot moves in indoor environment)
     *              set it to zero. Otherwise, for navigation on uneven terrains, it is normally modeled as 1DoF along the
     *              axis of the pitch angle of your robot.
     *
     * Note that it is very important how you specify the number of Trees according to the contact points.
     * It is assumed that one Tree can only have one single point of contact with the ground at the same time.
     * Therefore, if your real/physical kinematic Tree can have more than one contact point at a time (e.g: two)
     * the Tree needs to be split according to it (e.g: two Trees).
     * Take <a href="http://robotik.dfki-bremen.de/en/forschung/robotersysteme/asguard-ii.html">Asguard</a> wheel
     * as an example. If we want to model the wheel as two feet can have point in contact, two Trees
     * needs to be created in the wheel (virtually increasing the number of Trees).
     *
     */
    template <typename _Scalar>
    class MotionModel
    {
        protected:
            int number_trees, number_robot_joints, number_slip_joints, number_contact_joints;

        public:
            unsigned int model_dof;

        protected:

            /**@brief Forms the Navigation Equations for the navigation kinematics
             *
             * It assumes that known quantities are angular rotations
             * (Cartesian) and joint velocities.  non-slip in the X and Y axis.
             * Only slip in z-axis (non-holonomic constrain) The unknown
             * quantities are position of the robot in Cartesian, the contact
             * angle and the z value of the slip vector.  The system of
             * equations is unknownA * unknownx = knownB * knowny The
             * parameters are stored in the general dimension for the problem,
             * therefore unknown quantities are set to NaN in the parameters of
             * the method.
             *
             * @param[in] cartesianVelocities Eigen vector with the cartesian
             * velocities (w.r.t local body frame).
             *
             * @param[in] modelVelocities Eigen vector with the velocities of
             * the model(joint, slip and contact angle)
             *
             * @param[in] J robot jacobian matrix.
             *
             * @param[in] cartesianVelCov uncertainty in the measurement of the
             * cartesianVelocities vector.
             *
             * @param[in] modelVelCov uncertainty in the measurement of the
             * modelVelocities vector.
             *
             * @param[out] unknownA is the matrix of unknow quantities.
             * @param[out] unknownx in the vector of unknow quantities.
             * @param[out] knownB is the matrix of know quantities
             * @param[out] knowny is the matrix of know quantities
             * @param[out] Weight matrix with the inverse of the noise for the
             * equation (Weighted Least-Squares).
             */
            inline void navEquations(const Eigen::Matrix <_Scalar, 6, 1> &cartesianVelocities,
                                    const Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> &modelVelocities,
                                    const Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &J,
                                    const Eigen::Matrix <_Scalar, 6, 6> &cartesianVelCov,
                                    const Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &modelVelCov,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &unknownA,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> &unknownx,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &knownB,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> &knowny,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &Weight)
            {
                Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> spareI;
                spareI.resize(6*this->number_trees, 6);

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout<<"[MOTION_MODEL] navEquations\n";
                #endif

                /** Compute the composite rover equation matrix I **/
                for (register int i=0; i<this->number_trees; ++i)
                    spareI.block(i*6, 0, 6, 6) = Eigen::Matrix <_Scalar, 6, 6>::Identity();

                /** Form the unknownA matrix **/
                unknownA.block(0, 0, 6*this->number_trees, 3) = spareI.block(0, 0, 6*this->number_trees, 3);//Linear Velocities

                /** Form the unknownx vector **/
                unknownx.block(0, 0, 3, 1) = cartesianVelocities.block(0, 0, 3, 1);

                /** Non-holonomic constraint Slip vector **/
                unknownA.block(0, 3, 6*this->number_trees, this->number_trees) = -J.block(0, this->number_robot_joints+this->number_slip_joints-(this->number_trees), 6*this->number_trees, this->number_trees);
                unknownx.block(3, 0, this->number_trees, 1) = modelVelocities.block(this->number_robot_joints+this->number_slip_joints-(this->number_trees), 0, this->number_trees, 1);

                /** Contact angles per each Tree **/
                unknownA.block(0, 3+this->number_trees, 6*this->number_trees, this->number_contact_joints) = -J.block(0, this->number_robot_joints+this->number_slip_joints, 6*this->number_trees, this->number_contact_joints);
                unknownx.block(3+this->number_trees, 0, this->number_contact_joints, 1) = modelVelocities.block(this->number_robot_joints+this->number_slip_joints, 0, this->number_contact_joints, 1);

                /** Form the knownB matrix **/
                knownB.block(0, 0, 6* this->number_trees, 3) = -spareI.block(0, 3, 6*this->number_trees, 3); //Angular velocities
                knownB.block(0, 3, 6* this->number_trees, this->number_robot_joints) = J.block(0, 0, 6*this->number_trees, this->number_robot_joints);

                /** Form the knowny vector **/
                knowny.block(0, 0, 3, 1) = cartesianVelocities.template block<3, 1> (3,0);
                knowny.block(3, 0, this->number_robot_joints, 1) = modelVelocities.block(0, 0, this->number_robot_joints, 1);

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout<< "[MOTION_MODEL] spareI is of size "<<spareI.rows()<<"x"<<spareI.cols()<<"\n";
                std::cout<< "[MOTION_MODEL] The spareI matrix:\n" << spareI << std::endl;
                std::cout<< "[MOTION_MODEL] unknownA is of size "<<unknownA.rows()<<"x"<<unknownA.cols()<<"\n";
                std::cout<< "[MOTION_MODEL] The unknownA matrix:\n" << unknownA << std::endl;
                std::cout<< "[MOTION_MODEL] knownB is of size "<<knownB.rows()<<"x"<<knownB.cols()<<"\n";
                std::cout<< "[MOTION_MODEL] The knownB matrix:\n" << knownB << std::endl;
                std::cout<< "[MOTION_MODEL] Weight is of size "<<Weight.rows()<<"x"<<Weight.cols()<<"\n";
                std::cout<< "[MOTION_MODEL] The Weight matrix:\n" << Weight << std::endl;
                #endif

                return;
            }

        public:

            /**@brief Constructor
             *
             * Constructs the Motion Model object with the parameters.
             *
             * @param[out] status, true if everything went well.
             * @param[in] method of the points in contact selection algorithm
             * @param[in] pointer to the object with the Kinematic Model of the robot.
             *
             */
            MotionModel (const int _number_trees, const int _number_robot_joints, const int _number_slip_joints, const int _number_contact_joints):
                            number_trees(_number_trees),
                            number_robot_joints(_number_robot_joints),
                            number_slip_joints(_number_slip_joints),
                            number_contact_joints(_number_contact_joints)
            {
                this->model_dof = number_robot_joints + number_slip_joints + number_contact_joints;

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout<<"[MOTION_MODEL] Constructor. Model DoF:"<<this->model_dof<<"\n";
                #endif

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout<<"\n[MOTION_MODEL] **** END ****\n";
                #endif

                return;
            }


            /**@brief Default destructor
             */
            ~MotionModel()
            {
            }

            /**@bief Solver for the Navigation equations
             *
             * This method get the robot Jacobian matrix and computes the Navigation
             * kinematics equations to finally solve the Weighting Least-Square method.
             *
             * The Navigation kinematics equation: A * x = B * y
             * This method(navigation solver) perform the solution to the system above, solving:
             *
             *               x = (A^T * W * A)^(-1) * A^T * W * B * y
             *
             * @param[in] modelPositions, position values of the model (joints, slip and contact angle)
             * @param[in, out] cartesianVelocities, velocities of the robot Cartesian space (w.r.t local body frame).
             * @param[in] modelVelocities, velocity values of the model 9joints, slip and contact angle)
             * @param[in] cartesianVelCov, uncertainty in the measurement of the robot (local body frame) Cartesian velocities.
             * @param[in] modelVelCov, uncertainty in the measurement of the robot model velocities (joint space).
             * @param[in] Weight, trees Weighting matrix (which robot's tree has more weight to the computation).
             */
            virtual double navSolver(const Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> &modelPositions,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> &modelVelocities,
                                    const Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &J,
                                    Eigen::Matrix <_Scalar, 6, 1> &cartesianVelocities,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> &modelVelCov,
                                    Eigen::Matrix <_Scalar, 6, 6> &cartesianVelCov,
                                    Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> Weight)
            {
                double normalizedError = std::numeric_limits<double>::quiet_NaN(); //solution error of the Least-Squares
                std::vector<_Scalar> vectorPositions; // model positions in std_vector form
                vectorPositions.resize(model_dof);
                std::fill(vectorPositions.begin(), vectorPositions.end(), 0);
                Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> unknownA; // Nav non-sensed values matrix
                unknownA.resize(6*this->number_trees, 3+this->number_trees+this->number_contact_joints);
                Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> unknownx; // Nav non-sensed values vector
                unknownx.resize(3+this->number_trees+this->number_contact_joints, 1);
                Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> knownB; // Nav sensed values matrix
                knownB.resize(6*this->number_trees, 3+this->number_robot_joints);
                Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> knowny; // Nav sensed values vector
                knowny.resize(3+this->number_robot_joints, 1);

                /** Initialize  variables **/
                unknownA.setZero(); unknownx.setZero();
                knownB.setZero(); knowny.setZero();

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout << "[MOTION_MODEL] cartesianVelocities is of size "<<cartesianVelocities.rows()<<"x"<<cartesianVelocities.cols()<<"\n";
                std::cout << "[MOTION_MODEL] cartesianVelocities is \n" << cartesianVelocities<< std::endl;

                std::cout << "[MOTION_MODEL] modelVelocities is of size "<<modelVelocities.rows()<<"x"<<modelVelocities.cols()<<"\n";
                std::cout << "[MOTION_MODEL] modelVelocities is \n" << modelVelocities<< std::endl;
                #endif

                /** Assert vectors sizes **/
                assert(modelPositions.size() == vectorPositions.size());
                assert(modelPositions.size() == modelVelocities.size());
                assert(modelVelCov.cols() == modelVelCov.rows());
                assert(modelVelCov.cols() == modelVelocities.size());
                assert(Weight.cols() == 6*this->number_trees);
                assert(Weight.rows() == 6*this->number_trees);

                /** Copy Eigen to vector **/
                Eigen::Map <Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> > (&(vectorPositions[0]), model_dof) = modelPositions;

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout<< "[MOTION_MODEL] J is of size "<<J.rows()<<"x"<<J.cols()<<"\n";
                std::cout<< "[MOTION_MODEL] The J matrix \n" << J << std::endl;
                #endif

                /** Form the Composite Navigation Equations and Noise Covariance **/
                this->navEquations (cartesianVelocities, modelVelocities, J,
                        cartesianVelCov, modelVelCov, unknownA, unknownx, knownB, knowny, Weight);

                /** Solve the Motion Model by Least-Squares (navigation kinematics) **/
                Eigen::Matrix <_Scalar, Eigen::Dynamic, 1> knownb;
                knownb.resize(6*this->number_trees, 1);
                knownb = knownB*knowny;

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                /** DEBUG OUTPUT **/
                Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> Conj;
                Conj.resize(6*this->number_trees, 3+this->number_trees+this->number_contact_joints+1);

                Eigen::FullPivLU<Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> > lu_decompA(unknownA);
                std::cout << "[MOTION_MODEL] The rank of A is " << lu_decompA.rank() << std::endl;

                Eigen::FullPivLU<Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> > lu_decompB(knownB);
                std::cout << "[MOTION_MODEL] The rank of B is " << lu_decompB.rank() << std::endl;

                Conj.block(0, 0, 6*this->number_trees, 3+this->number_trees+this->number_contact_joints) = unknownA;
                Conj.block(0, 3+this->number_trees+this->number_contact_joints, 6*this->number_trees, 1) = knownb;
                Eigen::FullPivLU< Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> > lu_decompConj(Conj);
                std::cout << "[MOTION_MODEL] The rank of A|B*y is " << lu_decompConj.rank() << std::endl;
                std::cout << "[MOTION_MODEL] Pseudoinverse of A\n" << (unknownA.transpose() * Weight * unknownA).inverse() << std::endl;
                /*******************/
                #endif

                unknownx = (unknownA.transpose() * Weight * unknownA).ldlt().solve(unknownA.transpose() * Weight * knownb);

                /** Error of the solution **/
                Eigen::Matrix<double, 1,1> squaredError = (((unknownA*unknownx - knownb).transpose() * Weight * (unknownA*unknownx - knownb)));
                if (knownb.norm() != 0.00)
                    normalizedError = sqrt(squaredError[0]) / knownb.norm();


                /** Error covariance matrix **/
                Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> errorCov;
                errorCov.resize(6*this->number_trees, 6*this->number_trees);
                errorCov = (unknownA*unknownx - knownb).asDiagonal(); errorCov *= errorCov;// L-S error covariance
                Eigen::Matrix <_Scalar, Eigen::Dynamic, Eigen::Dynamic> uncertaintyCov; // noise cov
                uncertaintyCov.resize(3+this->number_trees+this->number_contact_joints, 3+this->number_trees+this->number_contact_joints);
                uncertaintyCov = (unknownA.transpose() * errorCov.inverse() * unknownA).inverse(); // Observer
                uncertaintyCov = 0.5*(uncertaintyCov + uncertaintyCov.transpose());// Guarantee symmetry

                /** Save the results in the parameters (previous NaN values are now just known quantities) **/
                cartesianVelocities.block(0, 0, 3, 1) = unknownx.block(0, 0, 3, 1); // Linear velocities
                cartesianVelCov.block(0, 0, 3, 3) = uncertaintyCov.block(0, 0, 3,3);//Linear Velocities noise

                /** Angular velocity estimated noise (experimental). Get the uncertainty of the angular velocity from the LS covariance **/
                if (cartesianVelCov.block(3, 3, 3, 3) == Eigen::Matrix3d::Zero())
                {
                    /** There is one per each _RobotTrees. errorCov is a 6*_RobotTrees x 6*_RobotTrees matrix dimension **/
                    for (register size_t i=0; i<this->number_trees; ++i)
                    {
                        cartesianVelCov.block(3, 3, 3, 3) += Weight.block(3+(6*i), 3+(6*i), 3, 3) * errorCov.block(3+(6*i),3+(6*i), 3, 3);//Angular Velocities noise
                    }
                }

                /** Non-holonomic constraint Slip vector **/
                modelVelocities.block(this->number_robot_joints+this->number_slip_joints-(this->number_trees), 0, this->number_trees, 1) = unknownx.block(3, 0, this->number_trees, 1);
                modelVelCov.block(this->number_robot_joints+this->number_slip_joints-(this->number_trees), this->number_robot_joints+this->number_slip_joints-(this->number_trees), this->number_trees, this->number_trees) = uncertaintyCov.block(3, 3, this->number_trees, this->number_trees);//pseudoInvUnknownA.col(3+i)[3+i]; For the time being set the error to the error in the estimation

                /** Contact angles per each Tree **/
                for (register int i=0; i<this->number_contact_joints; ++i)
                {
                    modelVelocities[this->number_robot_joints+this->number_slip_joints+i] = unknownx[3+this->number_trees+i];
                    modelVelCov.col(this->number_robot_joints+this->number_slip_joints+i)[this->number_robot_joints+this->number_slip_joints+i] = uncertaintyCov.col(3+this->number_trees+i)[3+this->number_trees+i];//pseudoInvUnknownA.col(3+_RobotTrees+i)[3+_RobotTrees+i];
                }

                #ifdef DEBUG_PRINTS_ODOMETRY_MOTION_MODEL
                std::cout << "[MOTION_MODEL] L-S solution:\n"<<unknownx<<std::endl;

                std::cout << "[MOTION_MODEL] RESULT errorCov is of size "<<errorCov.rows()<<"x"<<errorCov.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT errorCov is \n" << errorCov << std::endl;

                std::cout << "[MOTION_MODEL] RESULT uncertaintyCov is of size "<<uncertaintyCov.rows()<<"x"<<uncertaintyCov.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT uncertaintyCov is \n" << uncertaintyCov << std::endl;

                std::cout << "[MOTION_MODEL] RESULT cartesianVelocities is of size "<<cartesianVelocities.rows()<<"x"<<cartesianVelocities.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT cartesianVelocities is \n" << cartesianVelocities<< std::endl;

                std::cout << "[MOTION_MODEL] RESULT cartesianVelCov is of size "<<cartesianVelCov.rows()<<"x"<<cartesianVelCov.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT cartesianVelCov is \n" << cartesianVelCov<< std::endl;

                std::cout << "[MOTION_MODEL] RESULT modelVelocities is of size "<<modelVelocities.rows()<<"x"<<modelVelocities.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT modelVelocities is \n" << modelVelocities<< std::endl;

                std::cout << "[MOTION_MODEL] RESULT modelVelCov is of size "<<modelVelCov.rows()<<"x"<<modelVelCov.cols()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT modelVelCov is \n" << modelVelCov<< std::endl;

                std::cout << "[MOTION_MODEL] RESULT The absolute least squared error is:\n" << squaredError << std::endl;
                std::cout << "[MOTION_MODEL] RESULT The relative error is:\n" << normalizedError << std::endl;
                std::cout << "[MOTION_MODEL] RESULT The error vector is \n"<<(unknownA*unknownx - knownb)<<"\n";
                std::cout << "[MOTION_MODEL] RESULT The error variance is \n"<<errorCov<<"\n";
                std::cout << "[MOTION_MODEL] RESULT The error variance.inverse() is \n"<<errorCov.inverse()<<"\n";
                std::cout << "[MOTION_MODEL] RESULT The solution covariance is \n"<<cartesianVelCov<<"\n";
                #endif


                return normalizedError;
            }

             /**@brief Solver for the Slip equations
             * TO-DO
             */
            virtual base::Vector6d slipSolver(void)
            {
                base::Vector6d s;

                return s;
            }
    };
}

#endif //ODOMETRY_MOTION_MODEL_HPP

