/*
 * Copyright 2010,
 *
 * Francois Keith
 * Layale Saab
 * Olivier Stasse,
 *
 *
 * JRL/LAAS, CNRS/AIST
 *
 * This file is part of dynamicsJRLJapan.
 * dynamicsJRLJapan is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * dynamicsJRLJapan is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with dynamicsJRLJapan.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Research carried out within the scope of the Associated
 *  International Laboratory: Joint Japanese-French Robotics
 *  Laboratory (JRL)
 *
 */
#include <iostream>
#include <string>
#include <cstdio>
#include "jrl/dynamics/dynamicsfactory.hh"

using namespace std;

void DisplayMatrix4x4(MAL_S4x4_MATRIX(todisplay,double), ostream &os)
{
  for(unsigned int i=0;i<4;i++)
    {
      for(unsigned int j=0;j<4;j++)
	{
	  os << MAL_S4x4_MATRIX_ACCESS_I_J(todisplay,i,j) << " " ;
	}
      os << std::endl;
    }
}

void RecursiveDisplayOfJoints(CjrlJoint *aJoint, unsigned int verbosedisplay=0)
{
  if (aJoint==0)
    return;

  int NbChildren = aJoint->countChildJoints();

  std::cout << "CurrentTransformation " <<
    aJoint->currentTransformation() << std::endl;

  if (verbosedisplay>2)
    {
      std::cout << "Number of child  :" << NbChildren << std::endl;
      for(int i=0;i<NbChildren;i++)
	{
	  aJoint = aJoint->childJoint(i);
	}


      std::cout << "Nb of degree of freedom " <<
	aJoint->numberDof() << std::endl;

      MAL_S4x4_MATRIX_TYPE(double) initialpos = aJoint->initialPosition();
      std::cout << "Initial Position ";
      DisplayMatrix4x4(initialpos,std::cout);

      MAL_S4x4_MATRIX_TYPE(double) currentTransformation = aJoint->currentTransformation();
      std::cout << "CurrentTransformation ";
      DisplayMatrix4x4(currentTransformation,std::cout);

      std::cout << " Joint from root to here:" << std::endl;
      std::vector<CjrlJoint*> JointsFromRootToHere = aJoint->jointsFromRootToThis();

      std::cout << " Nb of nodes: " << JointsFromRootToHere.size() << std::endl;

      CjrlRigidVelocity aRV = aJoint->jointVelocity();
      std::cout << " Linear Velocity " << aRV.linearVelocity() << std::endl;
      std::cout << " Angular Velocity " << aRV.rotationVelocity() << std::endl;
      CjrlRigidAcceleration aRA = aJoint->jointAcceleration();
      std::cout << " Linear Acceleration " << aRA.linearAcceleration() << std::endl;
      std::cout << " Angular Acceleration " << aRA.rotationAcceleration() << std::endl;

      std::cout << "***********************************************" << std::endl;
      std::cout << " Display Now information related to children :" << std::endl;
    }

  for(int i=0;i<NbChildren;i++)
    {
      // Returns a const so we have to force the casting/
      RecursiveDisplayOfJoints((CjrlJoint *)aJoint->childJoint(i));
    }

}


void DisplayDynamicRobotInformation(CjrlDynamicRobot *aDynamicRobot)
{
  std::vector<CjrlJoint *> aVec = aDynamicRobot->jointVector();
  int r = aVec.size();
  std::cout << "Number of joints :" << r << std::endl;
}

void DisplayMatrix(MAL_MATRIX_TYPE(double) &aJ)
{
  for(int i=0;i<6;i++)
    {
      for(unsigned int j=0;j<MAL_MATRIX_NB_COLS(aJ);j++)
	{
	  if (aJ(i,j)==0.0)
	    printf("0 ");
	  else
	    printf("%10.5f ",aJ(i,j));
	}
      printf("\n");
    }

}

int main(int argc, char *argv[])
{

  if (argc!=6)
    {
      cerr << " This program takes 5 arguments: " << std::endl;
      cerr << "./TestHumanoidDynamicRobot PATH_TO_VRML_FILE VRML_FILE_NAME "<< std::endl;
      cerr << " PATH_TO_SPECIFICITIES_XML PATH PATH_TO_MAP_JOINT_2_RANK REDUCED_PATH_TO_MAP_JOINT_2_RANK" << std::endl;
      exit(-1);
    }


  string aSpecificitiesFileName = argv[3];
  string aPath=argv[1];
  string aName=argv[2];
  string aMapFromJointToRank=argv[4];
  string aMapFromJointToRankSmall=argv[5];

  CjrlJoint* rootJoint=0;

  dynamicsJRLJapan::ObjectFactory dynFactory;
  CjrlHumanoidDynamicRobot * aHDR  = dynFactory.createHumanoidDynamicRobot();
  string RobotFileName = aPath + aName;
  dynamicsJRLJapan::parseOpenHRPVRMLFile(*aHDR,RobotFileName,
					 aMapFromJointToRank,aSpecificitiesFileName);


  CjrlHumanoidDynamicRobot * aHDRSmall =  dynFactory.createHumanoidDynamicRobot();
  dynamicsJRLJapan::parseOpenHRPVRMLFile(*aHDRSmall,RobotFileName,
					 aMapFromJointToRank,aSpecificitiesFileName);

  // Display tree of the joints.
  std::cout << "Small model" << std::endl;
  rootJoint = aHDRSmall->rootJoint();

  // Test the tree.
  RecursiveDisplayOfJoints(rootJoint);


  // Tes the computation of the jacobian.
  double dInitPos[40] = {
    0.0, 0.0, -26.0, 50.0, -24.0, 0.0, 0.0, 0.0, -26.0, 50.0, -24.0, 0.0,  // legs

    0.0, 0.0, 0.0, 0.0, // chest and head

    15.0, -10.0, 0.0, -30.0, 0.0, 0.0, 10.0, // right arm
    15.0,  10.0, 0.0, -30.0, 0.0, 0.0, 10.0, // left arm

    -20.0, 20.0, -20.0, 20.0, -20.0, // right hand
    -10.0, 10.0, -10.0, 10.0, -10.0  // left hand
  };

  int NbOfDofs = aHDR->numberDof();
  std::cout << "NbOfDofs :" << NbOfDofs << std::endl;
  MAL_VECTOR_DIM(aCurrentConf,double,NbOfDofs);
  int lindex=0;
  for(int i=0;i<6;i++)
    aCurrentConf[lindex++] = 0.0;

  for(int i=0;i<(NbOfDofs-6 < 40 ? NbOfDofs-6 : 40) ;i++)
    aCurrentConf[lindex++] = dInitPos[i]*M_PI/180.0;
  //aCurrentConf[lindex++] = 0.0;

  aHDR->currentConfiguration(aCurrentConf);

  int NbOfDofsSmall = aHDRSmall->numberDof();
  std::cout << "NbOfDofs :" << NbOfDofsSmall << std::endl;
  MAL_VECTOR_DIM(aCurrentConfSmall,double,NbOfDofsSmall);
  lindex=0;
  for(int i=0;i<6;i++)
    aCurrentConfSmall[lindex++] = 0.0;

  for(int i=0;i<(NbOfDofsSmall-6 < 40 ? NbOfDofsSmall-6 : 40) ;i++)
    aCurrentConfSmall[lindex++] = dInitPos[i]*M_PI/180.0;
  //aCurrentConf[lindex++] = 0.0;

  aHDRSmall->currentConfiguration(aCurrentConfSmall);

  MAL_VECTOR_DIM(aCurrentVel,double,NbOfDofs);
  MAL_VECTOR_DIM(aCurrentVelSmall,double,NbOfDofsSmall);
  lindex=0;
  for(int i=0;i<NbOfDofsSmall;i++)
    aCurrentVelSmall[i] = aCurrentVel[i] =0.0;

  MAL_S3_VECTOR(ZMPval,double);

  aHDR->currentVelocity(aCurrentVel);
  aHDRSmall->currentVelocity(aCurrentVelSmall);

  //  aHDR->setComputeZMP(true);
  string inProperty="ComputeZMP"; string aValue="true";
  aHDR->setProperty(inProperty,aValue);
  aHDR->computeForwardKinematics();

  aHDRSmall->setProperty(inProperty,aValue);
  aHDRSmall->computeForwardKinematics();

  // Display tree of the joints.
  rootJoint = aHDR->rootJoint();

  // Test the tree.
  std::cout << "Normal model" << std::endl;
  RecursiveDisplayOfJoints(rootJoint);

  // Display tree of the joints.
  std::cout << "Small model" << std::endl;
  rootJoint = aHDRSmall->rootJoint();

  // Test the tree.
  RecursiveDisplayOfJoints(rootJoint);

  ZMPval = aHDR->zeroMomentumPoint();
  std::cout << "Normal model" << std::endl;
  std::cout << "Value of ZMP : " << ZMPval <<std::endl;
  std::cout << "Should be equal to the CoM: " << aHDR->positionCenterOfMass() << std::endl;

  ZMPval = aHDRSmall->zeroMomentumPoint();
  std::cout << "Small model" << std::endl;
  std::cout << "Value of ZMP : " << ZMPval <<std::endl;
  std::cout << "Should be equal to the CoM: " << aHDRSmall->positionCenterOfMass() << std::endl;

  std::vector<CjrlJoint *> aVec = aHDRSmall->jointVector();
  aVec[22]->computeJacobianJointWrtConfig();

  MAL_MATRIX_TYPE(double) aJ = aVec[22]->jacobianJointWrtConfig();

  DisplayMatrix(aJ);

  delete aHDR;
  delete aHDRSmall;

}
