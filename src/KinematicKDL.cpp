#include "KinematicKDL.hpp"

#ifndef D2R
#define D2R M_PI/180.00 /** Convert degree to radian **/
#endif
#ifndef R2D
#define R2D 180.00/M_PI /** Convert radian to degree **/
#endif

#define DEBUG_PRINTS 1

using namespace threed_odometry;

KinematicKDL::KinematicKDL (const int _RobotTrees, const int _RobotJointDoF,
                            const int _SlipDoF, const int _ContactDoF,
                            std::string &urdf_file, std::vector<float> &wheel_radius):
                            KinematicModel (_RobotTrees, _RobotJointDoF,  _SlipDoF, _ContactDoF)
{

    /** Take the wheel radius **/
    this->wheel_radius = wheel_radius;

    /** The number of wheels must be the same as the number of trees **/
    /** In case your robot does not have wheels, make a vector of wheels with zero radius **/
    assert (this->wheel_radius.size() == static_cast<size_t>(this->getNumberOfTrees()));

    /** Properly size the Jacobian matrices. One per wheel/One foot per wheel **/
    vector_jacobians.resize(this->wheel_radius.size());
    for(register size_t i=0; i<wheel_radius.size(); ++i)
    {
        vector_jacobians[i].resize(2*3, MAX_CHAIN_DOF);
    }

    /** Read URDF file **/
    std::string xml_string;
    const char * urdf_char = urdf_file.c_str();
    std::fstream xml_file(urdf_char, std::fstream::in);
    while ( xml_file.good() )
    {
        std::string line;
        std::getline( xml_file, line);
        xml_string += (line + "\n");
    }
    xml_file.close();

    boost::shared_ptr<urdf::ModelInterface> robot = urdf::parseURDF(xml_string);
    if (!robot)
    {
        throw std::runtime_error("[KDL_MODEL] Constructor could not parse URDF model\n");
    }

    /** Assign an information string name **/
    this->model_name = robot->getName();

    LOG_INFO("[KDL_MODEL] Robot name is: %s\n", this->model_name.c_str());
    LOG_INFO("[KDL_MODEL] Robot has %d number of Trees\n",this->wheel_radius.size());

    /* get root link*/
    boost::shared_ptr<const urdf::Link> root_link = robot->getRoot();
    if (!root_link)
        throw std::runtime_error("[KDL_MODEL] Constructor could not find Root link\n");

    LOG_INFO("[KDL_MODEL] Robot Link: %s has %d child(ren)\n", root_link->name.c_str(), root_link->child_links.size());

    if (!kdl_parser::treeFromUrdfModel(*robot, this->tree))
        throw std::runtime_error("[KDL_MODEL] Constructor could not extract KDL tree\n");

    KDL::SegmentMap::const_iterator root = this->tree.getRootSegment();
    #ifdef DEBUG_PRINTS
    std::cout << " ======================================" << std::endl;
    std::cout << " Original Tree has " << this->tree.getNrOfSegments() << " link(s) and "<< this->tree.getNrOfJoints() <<" Joint(s)\n";
    std::cout << " ======================================" << std::endl;
    #endif

    for (register unsigned int j=0; j<this->wheel_radius.size(); ++j)
    {
        std::stringstream chainidx;
        chainidx << j;
        #ifdef DEBUG_PRINTS
        std::cout<<"[KDL_MODEL] chainname: "<< chainidx.str()<<"\n";
        #endif

        KDL::Chain chain;//!Create the chain

        /** Create the contact angle and the slip vector/displacement subchain **/
        chain.addSegment(KDL::Segment(chainidx.str()+"_CONTACT", KDL::Joint(KDL::Joint::RotY),KDL::Frame(KDL::Vector(0.0,0.0,-wheel_radius[j]))));//Contact angle
        chain.addSegment(KDL::Segment(chainidx.str()+"_MOVE",KDL::Joint(KDL::Joint::TransX),KDL::Frame(KDL::Vector(0.0,0.0,0.0))));//Translation along X-axis
        chain.addSegment(KDL::Segment(chainidx.str()+"_SLIPX",KDL::Joint(KDL::Joint::TransX),KDL::Frame(KDL::Vector(0.0,0.0,0.0))));//Slip vector along X-axis
        chain.addSegment(KDL::Segment(chainidx.str()+"_SLIPY",KDL::Joint(KDL::Joint::TransY),KDL::Frame(KDL::Vector(0.0,0.0,0.0))));//Slip vector along Y-axis
        chain.addSegment(KDL::Segment(chainidx.str()+"_SLIPZ",KDL::Joint(KDL::Joint::RotZ),KDL::Frame(KDL::Vector(0.0,0.0,0.0))));//Slip vector along Z-axis

        /** Add to the end of the selected chain **/
        std::string endChain;
        switch (j)
        {
            case 1:
                 endChain = getContactPoint(root->second.children[1], 0);
                 break;
            case 2:
                 endChain = getContactPoint(root->second.children[0], 1);
                 break;
            case 3:
                 endChain = getContactPoint(root->second.children[1], 1);
                 break;
            case 4:
                 endChain = getContactPoint(root->second.children[2], 0);
                 break;
            case 5:
                 endChain = getContactPoint(root->second.children[2], 1);
                 break;
            default:
                 endChain = getContactPoint(root->second.children[0], 0);
                 break;

        }
        bool exit_value = tree.addChain(chain, endChain);
        if(!exit_value)
        {
            LOG_INFO("[KDL_MODEL] ERROR in Constructor\n");
        }
        LOG_INFO("[KDL_MODEL] Odometry Kinematics Chain has: %d new more segments\n", chain.getNrOfSegments());


    }

    #ifdef DEBUG_PRINTS
    /** Walk through tree **/
    std::cout << " ======================================" << std::endl;
    std::cout << " Tree has " << this->tree.getNrOfSegments() << " link(s) and "<< this->tree.getNrOfJoints() <<" Joint(s)\n";
    std::cout << " Root link" << std::endl;
    std::cout << " ======================================" << std::endl;
    printLink(root, "");
    #endif
}

KinematicKDL::~KinematicKDL()
{
}

std::string KinematicKDL::name ()
{
        return this->model_name;
}

void KinematicKDL::contactPointsPerTree (std::vector<unsigned int> &contactPoints)
{
    /** TO-DO **/
    /**  The robot has the same number of contact points per tree **/
    contactPoints.resize(this->wheel_radius.size());
    for (std::vector<unsigned int>::iterator it = contactPoints.begin() ; it != contactPoints.end(); ++it)
    {
        *it = 1;
    }

    return;
}

void KinematicKDL::setPointsInContact (const std::vector<int> &pointsInContact)
{
        this->contact_idx = pointsInContact;
}

void KinematicKDL::fkBody2ContactPointt(const int chainIdx, const std::vector<double> &positions,  Eigen::Affine3d &fkTrans, base::Matrix6d &fkCov)
{
    /** Check if the number of values is correct (Since a priori we don't know
     * which kinematics chain is, the maxchainDOF is needed) **/
    if (positions.size() == this->getMaxChainDoF())
    {
        KDL::SegmentMap::const_iterator root = this->tree.getRootSegment();

        if (chainIdx > this->RobotTrees)
        {
            fkTrans.matrix() = std::numeric_limits<double>::quiet_NaN() * Eigen::Matrix<double, 4, 4>::Identity();

            /** No uncertainty provided **/
            fkCov = std::numeric_limits<double>::quiet_NaN() * base::Matrix6d::Identity();

            return;
        }
        else
        {
            std::string contactPoint;
            switch (chainIdx)
            {
                case 1:
                     contactPoint = getContactPoint(root->second.children[1], 0);
                     break;
                case 2:
                     contactPoint = getContactPoint(root->second.children[0], 1);
                     break;
                case 3:
                     contactPoint = getContactPoint(root->second.children[1], 1);
                     break;
                case 4:
                     contactPoint = getContactPoint(root->second.children[2], 0);
                     break;
                case 5:
                     contactPoint = getContactPoint(root->second.children[2], 1);
                     break;
                default:
                     contactPoint = getContactPoint(root->second.children[0], 0);
                     break;


            }

            #ifdef DEBUG_PRINTS
            std::cout<<"[FORWARD_KINEMATICS] ContactPoint is :"<<contactPoint<<"\n";
            #endif

            /** Get the kinematics chain from the body to the contact point **/
            KDL::Chain chain;
            bool exit_value = tree.getChain(root->second.segment.getName(), contactPoint, chain);

            if (!exit_value)
            {
                #ifdef DEBUG_PRINTS
                std::cout<<"[FORWARD_KINEMATICS] Chain no found\n";
                #endif
            }

            /** Create solver based on kinematic chain **/
            KDL::ChainFkSolverPos_recursive fksolver = KDL::ChainFkSolverPos_recursive(chain);

            /** Create joint array **/
            unsigned int nj = chain.getNrOfJoints();
            KDL::JntArray jointpositions = KDL::JntArray(nj);
            #ifdef DEBUG_PRINTS
            std::cout<<"[FORWARD_KINEMATICS] Number of Joints is: "<<nj<<"\n";
            #endif

            /** Fill the joints positions (TO-DO: improve the coding stile with a fork) **/
            switch (chainIdx)
            {
                case 1:
                    jointpositions(0) = positions[1];//passive
                    jointpositions(1) = -positions[1];//mimic
                    jointpositions(2) = positions[4];//wheel walking
                    jointpositions(3) = positions[10];//steer
                    jointpositions(4) = positions[positions.size()-1];//contact angle
                    jointpositions(5) = this->wheel_radius[1]*positions[14];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(6) = positions[this->RobotJointDoF];
                    jointpositions(7) = positions[this->RobotJointDoF+1];
                    jointpositions(8) = positions[this->RobotJointDoF+2];

                    break;
                case 2:
                    jointpositions(0) = positions[0];//passive
                    jointpositions(1) = -positions[0];//mimic
                    jointpositions(2) = positions[5];//wheel walking
                    jointpositions(3) = positions[positions.size()-1];//contact angle
                    jointpositions(4) = this->wheel_radius[2]*positions[15];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(5) = positions[this->RobotJointDoF];
                    jointpositions(6) = positions[this->RobotJointDoF+1];
                    jointpositions(7) = positions[this->RobotJointDoF+2];

                    break;
                case 3:
                    jointpositions(0) = positions[1];//passive
                    jointpositions(1) = -positions[1];//mimic
                    jointpositions(2) = positions[6];//wheel walking
                    jointpositions(3) = positions[positions.size()-1];//contact angle
                    jointpositions(4) = this->wheel_radius[3]*positions[16];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(5) = positions[this->RobotJointDoF];
                    jointpositions(6) = positions[this->RobotJointDoF+1];
                    jointpositions(7) = positions[this->RobotJointDoF+2];


                    break;
                case 4:
                    jointpositions(0) = positions[2];//passive
                    jointpositions(1) = -positions[2];//mimic
                    jointpositions(2) = positions[7];//wheel walking
                    jointpositions(3) = positions[11];//steer
                    jointpositions(4) = positions[positions.size()-1];//contact angle
                    jointpositions(5) = this->wheel_radius[4]*positions[17];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(6) = positions[this->RobotJointDoF];
                    jointpositions(7) = positions[this->RobotJointDoF+1];
                    jointpositions(8) = positions[this->RobotJointDoF+2];

                    break;
                case 5:
                    jointpositions(0) = positions[2];//passive
                    jointpositions(1) = -positions[2];//mimic
                    jointpositions(2) = positions[8];//wheel walking
                    jointpositions(3) = positions[12];//steer
                    jointpositions(4) = positions[positions.size()-1];//contact angle
                    jointpositions(5) = this->wheel_radius[5]*positions[18];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(6) = positions[this->RobotJointDoF];
                    jointpositions(7) = positions[this->RobotJointDoF+1];
                    jointpositions(8) = positions[this->RobotJointDoF+2];

                    break;
                default:
                    jointpositions(0) = positions[0];//passive
                    jointpositions(1) = -positions[0];//mimic
                    jointpositions(2) = positions[3];//wheel walking
                    jointpositions(3) = positions[9];//steer
                    jointpositions(4) = positions[positions.size()-1];//contact angle
                    jointpositions(5) = this->wheel_radius[0]*positions[13];//displacement due to wheel rotation

                    /** Slip vector **/
                    jointpositions(6) = positions[this->RobotJointDoF];
                    jointpositions(7) = positions[this->RobotJointDoF+1];
                    jointpositions(8) = positions[this->RobotJointDoF+2];

                    break;
            }


            /** Create the frame that will contain the results **/
            KDL::Frame cartpos;

            /** Calculate forward position kinematics **/
            bool kinematics_status;
            kinematics_status = fksolver.JntToCart(jointpositions,cartpos);
            if(kinematics_status>=0)
            {
                /** Store the solution in the arguments **/
                double x,y,z,w;
                cartpos.M.GetQuaternion(x,y,z,w);
                Eigen::Quaternion<double> orientation(w,x,y,z);
                fkTrans = Eigen::Affine3d(orientation);
                fkTrans.translation() = Eigen::Vector3d (cartpos.p.x(), cartpos.p.y(), cartpos.p.z());

            }

            /** No uncertainty provided **/
            fkCov.setZero();
        }
    }

    return;
}

void KinematicKDL::fkSolver(const std::vector<double> &positions, std::vector<Eigen::Affine3d> &fkRobot, std::vector<base::Matrix6d> &fkCov)
{
    register unsigned int chainit = 0;
    std::vector<double> chainpositions;

    /** Check if the number of values is correct **/
    if (positions.size() == KinematicKDL::MODEL_DOF)
    {
        /** Resize the vectors **/
        fkRobot.resize(this->wheel_radius.size());
        fkCov.resize(this->wheel_radius.size());

        /** Resize the chain positions **/
        chainpositions.resize(this->getMaxChainDoF());
        std::fill(chainpositions.begin(), chainpositions.end(), 0.00);

        /** Fill the beginning (common part) **/
        for (register int i=0; i<this->RobotJointDoF; ++i)
        {
            chainpositions[i] = positions[i];
        }

        /** Calculate the Forward Kinematics for all the chains of the robot **/
        for (register int i=0; i<static_cast<int>(this->wheel_radius.size()); ++i)
        {
            /** Create the joints for this chain (the rest) **/
            for (register int l=0; l<this->SlipDoF; ++l)
            {
                chainpositions[this->RobotJointDoF+l] = positions[this->RobotJointDoF+(this->SlipDoF*i)+l];
            }

            /** Contact angle **/
            chainpositions[chainpositions.size()-1] = positions[KinematicKDL::MODEL_DOF-(wheel_radius.size()-i)];

            #ifdef DEBUG_PRINTS
            std::cout << "[ROBOT_FKSOLVER] chainpositions contains:";
            for (std::vector<double>::iterator it = chainpositions.begin() ; it != chainpositions.end(); ++it)
                std::cout << ' ' << *it;
            std::cout << '\n';
            #endif

            /** Perform the forward kinematics **/
            fkBody2ContactPointt(chainit, chainpositions, fkRobot[chainit], fkCov[chainit]);

            chainit++;
        }
    }
    return;

}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> KinematicKDL::jacobianSolver(const std::vector<double> &positions)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J;
    J.resize(6*this->RobotTrees, this->RobotJointDoF);

    /** Check if the number of values is correct **/
    if (positions.size() == static_cast<unsigned int>(this->RobotJointDoF))
    {
        /** Calculate the Jacobian Solver for all the chains of the robot **/
        for (register int i=0; i<static_cast<int>(this->wheel_radius.size()); ++i)
        {
            /** Get the chain **/
            std::stringstream chainidx;
            chainidx << i;

            KDL::SegmentMap::const_iterator root = this->tree.getRootSegment();
            std::string rootName = root->second.segment.getName();
            KDL::Chain chain, chainCitoCit, chainToContactPoint, chainFromContactPoint;
            bool exit_value1, exit_value2, exit_value3;
            exit_value1 = tree.getChain(rootName, chainidx.str()+"_CONTACT", chainToContactPoint);
            exit_value2 = tree.getChain(chainidx.str()+"_CONTACT", chainidx.str()+"_SLIPZ", chainCitoCit); //!CONTACT link is not taken (KDL works like that!!!)
            exit_value3 = tree.getChain(chainidx.str()+"_CONTACT", rootName, chainFromContactPoint);
            if (!exit_value1 || !exit_value2 || !exit_value3)
            {
                #ifdef DEBUG_PRINTS
                std::cout<<"[JACOBIAN_SOLVER_KDL] ERROR Chain no found\n";
                #endif
            }
            else
            {
                #ifdef DEBUG_PRINTS
                std::cout<<"[JACOBIAN_SOLVER_KDL] ChainIdx from "<<chainidx.str()+"_SLIPZ to "<< rootName <<"\n";
                #endif
            }

            /** Connect the kinematics chains and the root Segment **/
            chain.addSegment(root->second.segment); //! base_link
            chain.addChain(chainToContactPoint); //! From Body To contact plane
            chain.addChain(chainCitoCit); //! Movement and slip vector
            chain.addChain(chainFromContactPoint); //! From Contact To Body base

            #ifdef DEBUG_PRINTS
            std::cout << " ======================================" << std::endl;
            std::cout << " Chain has " << chain.getNrOfSegments() << " number of segments\n";
            for (register int l = 0; l<static_cast<int>(chain.getNrOfSegments()); ++l)
                std::cout <<"Segment["<<l<<"]: "<<chain.getSegment(l).getName() <<"\n";
            std::cout << " ======================================" << std::endl;
            #endif

            /** Create joint array **/
            unsigned int njTotal = chain.getNrOfJoints();
            KDL::JntArray jointpositions = KDL::JntArray(njTotal);

            /** Create KDL jacobian **/
            KDL::Jacobian kdljacob;
            kdljacob.resize(njTotal);

            /** Fill the joints positions (TO-DO: improve the coding stile with a fork) **/
            switch (i)
            {
                case 1:
                    jointpositions(0) = positions[1];//passive
                    jointpositions(1) = -positions[1];//mimic
                    jointpositions(2) = positions[4];//wheel walking
                    jointpositions(3) = positions[10];//steer
                    jointpositions(4) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+1];
                    jointpositions(5) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(6+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;
                case 2:
                    jointpositions(0) = positions[0];//passive
                    jointpositions(1) = -positions[0];//mimic
                    jointpositions(2) = positions[5];//wheel walking
                    jointpositions(3) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+2];
                    jointpositions(4) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(5+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;
                case 3:
                    jointpositions(0) = positions[1];//passive
                    jointpositions(1) = -positions[1];//mimic
                    jointpositions(2) = positions[6];//wheel walking
                    jointpositions(3) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+3];
                    jointpositions(4) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(5+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;
                case 4:
                    jointpositions(0) = positions[2];//passive
                    jointpositions(1) = -positions[2];//mimic
                    jointpositions(2) = positions[7];//wheel walking
                    jointpositions(3) = positions[11];//steer
                    jointpositions(4) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+4];
                    jointpositions(5) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(6+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;
                case 5:
                    jointpositions(0) = positions[2];//passive
                    jointpositions(1) = -positions[2];//mimic
                    jointpositions(2) = positions[8];//wheel walking
                    jointpositions(3) = positions[12];//steer
                    jointpositions(4) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+5];
                    jointpositions(5) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(6+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;

                default:
                    jointpositions(0) = positions[0];//passive
                    jointpositions(1) = -positions[0];//mimic
                    jointpositions(2) = positions[3];//wheel walking
                    jointpositions(3) = positions[9];//steer
                    jointpositions(4) = positions[this->RobotJointDoF+(this->RobotTrees*this->SlipDoF)+0];
                    jointpositions(5) = 0.00;//displacement due to wheel rotation

                    /** Slip vector **/
                    for(register int j=0; j<static_cast<int>(this->SlipDoF); ++j)
                        jointpositions(6+j) = positions[this->RobotJointDoF+(i*this->SlipDoF)+j];

                    break;
            }

            /** The sub-chain from Ci to Body **/
            double auxSubChainNrJoint = chainToContactPoint.getNrOfJoints() + chainCitoCit.getNrOfJoints();
            for (register int j=0; j<static_cast<int>(chainFromContactPoint.getNrOfJoints()); ++j)
                jointpositions(auxSubChainNrJoint+j) = jointpositions(chainToContactPoint.getNrOfJoints()-1-j);

            #ifdef DEBUG_PRINTS
            std::cout<<"[JACOBIAN_SOLVER_KDL] njTotal is: "<<njTotal<<"\n";
            std::cout<<"[JACOBIAN_SOLVER_KDL] Joints to Ci: "<<chainToContactPoint.getNrOfJoints()<<"\n";
            std::cout<<"[JACOBIAN_SOLVER_KDL] Joints displacement: "<<chainCitoCit.getNrOfJoints()<<"\n";
            std::cout<<"[JACOBIAN_SOLVER_KDL] Joints from Ci: "<<chainFromContactPoint.getNrOfJoints()<<"\n";
            std::cout << "[JACOBIAN_SOLVER_KDL] jointpositions contains:";
            for (register unsigned int it = 0; it<njTotal; ++it)
                std::cout << ' ' << jointpositions(it);
            std::cout << '\n';
            #endif

            /** Solver for the Jacobian **/
            KDL::ChainJntToJacSolver jacsolver = KDL::ChainJntToJacSolver(chain);
            jacsolver.JntToJac(jointpositions, kdljacob);

            /** Get the Jacobian in Eigen Matrix class **/
            Eigen::Matrix<double, 6, Eigen::Dynamic> eigenjacob;
            eigenjacob.resize(6,njTotal);
            eigenjacob = kdljacob.data;
            #ifdef DEBUG_PRINTS
            std::cout<<"[JACOBIAN_SOLVER_KDL] Jacobian is "<<kdljacob.rows()<<" x " <<kdljacob.columns() <<std::endl;
            std::cout<<"[JACOBIAN_SOLVER_KDL] Jacobian:\n"<<eigenjacob<<std::endl;
            #endif

            /** Get the Jacobian in the desired form **/
            vector_jacobians[i].setZero();
            switch(i)
            {
                /** FR Wheel**/
             case 1:
                vector_jacobians[i].col(1) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(4) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(10) = eigenjacob.col(njTotal-4);//steer joint
                vector_jacobians[i].col(14) = this->wheel_radius[1] * eigenjacob.col(5);//wheel rotation
                break;

                /** ML Wheel**/
             case 2:
                vector_jacobians[i].col(0) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(5) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(15) = this->wheel_radius[2] * eigenjacob.col(4);//wheel rotation
                break;

                /** MR Wheel**/
             case 3:
                vector_jacobians[i].col(1) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(6) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(16) = this->wheel_radius[3] * eigenjacob.col(4);//wheel rotation
                break;

                /** RL Wheel**/
             case 4:
                vector_jacobians[i].col(2) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(7) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(11) = eigenjacob.col(njTotal-4);//steer joint
                vector_jacobians[i].col(17) = this->wheel_radius[4] * eigenjacob.col(5);//wheel rotation
                break;

                /** RR Wheel**/
             case 5:
                vector_jacobians[i].col(2) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(8) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(12) = eigenjacob.col(njTotal-4);//steer joint
                vector_jacobians[i].col(18) = this->wheel_radius[5] * eigenjacob.col(5);//wheel rotation
                break;

                /** FL Wheel**/
             default:
                vector_jacobians[i].col(0) = eigenjacob.col(njTotal-2) - eigenjacob.col(njTotal-1);//passive parallel joints
                vector_jacobians[i].col(3) = eigenjacob.col(njTotal-3);//wheel walking
                vector_jacobians[i].col(9) = eigenjacob.col(njTotal-4);//steer joint
                vector_jacobians[i].col(13) = this->wheel_radius[0] * eigenjacob.col(5);//wheel rotation
                break;
            }

            switch(i)
            {
             case 0: case 1: case 4: case 5:
                /** Slip Vector and Contact angles **/
                vector_jacobians[i].block(0, this->RobotJointDoF, 6, this->SlipDoF+this->ContactDoF) =
                    eigenjacob.block(0, 6, 6, this->SlipDoF+this->ContactDoF);
                break;
             case 2: case 3:
                /** Slip Vector and Contact angles **/
                vector_jacobians[i].block(0, this->RobotJointDoF, 6, this->SlipDoF+this->ContactDoF) =
                    eigenjacob.block(0, 5, 6, this->SlipDoF+this->ContactDoF);
                break;
             default:
                break;
            }


            #ifdef DEBUG_PRINTS
            std::cout<<"[JACOBIAN_SOLVER_KDL] vector_jacobians is "<<vector_jacobians[i].rows()<<" x "<<vector_jacobians[i].cols()<<std::endl;
            std::cout<<"[JACOBIAN_SOLVER_KDL] vector_jacobians:\n"<<vector_jacobians[i]<<std::endl;
            #endif

        }

        J.setZero();

        /** Form the robot Jacobian matrix (sparse Jacobian matrix) **/
        for (register int i=0; i<static_cast<int>(this->wheel_radius.size()); ++i)
        {
            /** Joints including wheel rotation **/
            J.block(6*i, 0, 6, this->RobotJointDoF) = vector_jacobians[i].block(0,0, 6, this->RobotJointDoF);
            /** Slip vectors **/
            J.block(6*i, this->RobotJointDoF+(this->SlipDoF*i), 6, this->SlipDoF) = vector_jacobians[i].block(0,this->RobotJointDoF, 6, this->SlipDoF);
            /** Contact Angles **/
            J.block(6*i, this->RobotJointDoF+(this->SlipDoF*this->RobotTrees)+(this->ContactDoF*i), 6, this->ContactDoF) = vector_jacobians[i].block(0,this->RobotJointDoF+this->SlipDoF,6, this->ContactDoF);
        }

    }
    return J;
}
