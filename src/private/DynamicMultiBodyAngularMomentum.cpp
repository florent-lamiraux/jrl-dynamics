/* @doc Computation of the dynamic aspect for a robot.
   This class will load the description of a robot from a VRML file
   following the OpenHRP syntax. Using ForwardVelocity it is then
   possible specifying the angular velocity and the angular value 
   to get the absolute position, and absolute velocity of each 
   body separetly. Heavy rewriting from the original source
   of Adrien and Jean-Remy. 
 
   This implantation is an updated based on a mixture between 
   the code provided by Jean-Remy and Adrien.
 
   Copyright (c) 2005-2006, 
   @author Olivier Stasse, Ramzi Sellouati, Jean-Remy Chardonnet, Adrien Escande, Abderrahmane Kheddar
   Copyright (c) 2007-2009
   @author Olivier Stasse, Oussama Kannoun, Fumio Kanehiro.
   JRL-Japan, CNRS/AIST
 
   All rights reserved.

   Please refers to file License.txt for details on the license.

*/

/*! System includes */
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include "Debug.h"

/*! Local library includes. */
#include "MatrixAbstractLayer/MatrixAbstractLayer.h"
#include "dynamicsJRLJapan/DynamicBody.h"
#include "DynMultiBodyPrivate.h"
#include "robotDynamics/jrlBody.h"

#include "fileReader.h"

using namespace dynamicsJRLJapan;
/*! Implements the angular Momentum methods of DynMultiBodyPrivate */

void DynMultiBodyPrivate::angularMomentumWrtCoM(vector3d & angularmomentum) 
{
  angularMomentumWrtToPt(positionCoMPondere, angularmomentum);
}

void DynMultiBodyPrivate::angularMomentumWrtToPt(vector3d &apoint, vector3d & angularmomentum)
{
  /** Intermediate variables. The mantra is :
      "To optimize those variables, in the Compiler we trust"
      (with the appropriate compilation options).
  */
  vector3d NE_lP,NE_lw_c, NE_tmp3, NE_tmp2, NE_tmp,NE_lL;
  matrix3d NE_Rtmp, NE_Rt, NE_Ro, NE_Rot;
  /* End of intermediate */

  DynamicBodyPrivate *aDB=0;
  int currentNode = labelTheRoot;
  currentNode = m_listOfBodies[labelTheRoot]->child;
  vector3d lL(0.0,0.0,0.0);

  do
    {

      aDB = m_listOfBodies[currentNode];

      NE_lP = m_listOfBodies[currentNode]->P;
      ODEBUG("P: " << NE_lP );
      NE_lw_c = m_listOfBodies[currentNode]->w_c - positionCoMPondere;

      // Computes angular momentum matrix L
      // Lk = xc x Pk + R * I * Rt * w
      MAL_S3x3_TRANSPOSE_A_in_At(m_listOfBodies[currentNode]->R,NE_Rt);

      MAL_S3_VECTOR_CROSS_PRODUCT(NE_tmp3,NE_lw_c,NE_lP);

      MAL_S3x3_C_eq_A_by_B(NE_tmp2,NE_Rt , m_listOfBodies[currentNode]->w);
      MAL_S3x3_C_eq_A_by_B(NE_tmp, m_listOfBodies[currentNode]->getInertie(),NE_tmp2);
      MAL_S3x3_C_eq_A_by_B(NE_tmp2, m_listOfBodies[currentNode]->R,NE_tmp);
      NE_lL = NE_tmp3 + NE_tmp2;
      ODEBUG("L: " << lL);

      lL += NE_lL;

      int step=0;
      int NextNode=0;
      do
        {

	  if (step==0)
            {
	      NextNode = m_listOfBodies[currentNode]->child;
	      step++;
            }
	  else if(step==1)
            {
	      NextNode = m_listOfBodies[currentNode]->sister;
	      step++;
            }
	  else if (step==2)
            {
	      NextNode = m_listOfBodies[currentNode]->getLabelMother();
	      if (NextNode>=0)
                {
		  /* Test if current node is leaf,
		     because in this case the force are not set properly. */
		  if (m_ComputeBackwardDynamics)
                    {
		      if ((m_listOfBodies[currentNode]->sister==-1) &&
			  (m_listOfBodies[currentNode]->child==-1))
			BackwardDynamics(*m_listOfBodies[currentNode]);

		      /* Compute backward dynamics */
		      BackwardDynamics(*m_listOfBodies[NextNode]);
                    }
		  currentNode = NextNode;
		  NextNode = m_listOfBodies[currentNode]->sister;
                }
	      else
		NextNode=labelTheRoot;
            }


        }
      while (NextNode==-1);
      currentNode = NextNode;

    }
  while(currentNode!=labelTheRoot);
  
  angularmomentum = lL;
}

/**
   \brief Get the angular momentum of the robot at the center of mass.
*/
const MAL_S3_VECTOR(,double)& DynMultiBodyPrivate::angularMomentumRobot()
{
  return m_L;

} ;

/**
   \brief Get the time-derivative of the angular momentum at the center of mass.
*/
const MAL_S3_VECTOR(,double)& DynMultiBodyPrivate::derivativeAngularMomentum()
{

  return m_dL;

};

MAL_S3_VECTOR(,double) DynMultiBodyPrivate::GetL(int JointID)
{
  MAL_S3_VECTOR(empty,double);
  if ((JointID>=0) &&
      ((unsigned int)JointID<m_listOfBodies.size()))
    return m_listOfBodies[ConvertIDInActuatedToBodyID[JointID]]->L;
  return empty;
}

void DynMultiBodyPrivate::getJacobianAngularMomentumWrtCoM(matrixNxP &outjacobian)
{
  if ((MAL_MATRIX_NB_ROWS(outjacobian) != 3) || 
      (MAL_MATRIX_NB_COLS(outjacobian) != numberDof()))
    MAL_MATRIX_RESIZE(outjacobian,3,numberDof());

  MAL_MATRIX_FILL(outjacobian,0);

  unsigned int rank;
  JointPrivate* aJoint;
  DynamicBodyPrivate* aBody;
    
  for(unsigned int i=0;i<m_ConfigurationToJoints.size();i++)
    {
      if (m_ConfigurationToJoints[i] == rootJoint())
	continue;

      aJoint = m_ConfigurationToJoints[i];
      aBody=  aJoint->linkedDBody();
      rank = aJoint->rankInConfiguration();
      
      matrixNxP pJacobian;
      vector3d av(0,0,0); // Dummy 
      MAL_MATRIX_RESIZE(pJacobian,6, numberDof());
      getJacobian(*rootJoint(),*aJoint,av,pJacobian,true);

      ODEBUG("pJacobian:" <<pJacobian);
      matrixNxP pLinearJacobian;
      MAL_MATRIX_RESIZE(pLinearJacobian,3,MAL_MATRIX_NB_COLS(pJacobian));
      MAL_MATRIX_C_eq_EXTRACT_A(pLinearJacobian,pJacobian,double,0,0,3,
				MAL_MATRIX_NB_COLS(pJacobian));
      ODEBUG("pLinearJacobian:" <<endl <<pLinearJacobian);

      matrixNxP pAngularJacobian; 
      MAL_MATRIX_RESIZE(pAngularJacobian,3,MAL_MATRIX_NB_COLS(pJacobian));
      MAL_MATRIX_C_eq_EXTRACT_A(pAngularJacobian,pJacobian,double,3,0,3,
				MAL_MATRIX_NB_COLS(pJacobian));

      ODEBUG("pAngularJacobian:" <<endl <<pAngularJacobian);

      // Used to compute the anti-symmetric matrix.
      matrixNxP xkmxg_cp;double lmasse = aBody->getMasse();
      MAL_MATRIX_RESIZE(xkmxg_cp,3,3);
      av =aBody->w_c - positionCoMPondere;
      xkmxg_cp(0,0) =          0.0; xkmxg_cp(0,1) = -lmasse*av(2); xkmxg_cp(0,2) = lmasse*av(1);
      xkmxg_cp(1,0) = lmasse*av(2); xkmxg_cp(1,1) =           0.0; xkmxg_cp(1,2) =-lmasse*av(0);
      xkmxg_cp(2,0) =-lmasse*av(1); xkmxg_cp(2,1) =  lmasse*av(0); xkmxg_cp(2,2) =         0.0;

      ODEBUG("xkmxg_cp: " <<xkmxg_cp);

      matrixNxP leftoperand;
      MAL_C_eq_A_by_B(leftoperand,xkmxg_cp,pLinearJacobian);
      outjacobian = outjacobian + leftoperand;
      
      matrixNxP rightoperand;
      matrix3d tmp2_3d;
      matrixNxP tmp2;
      MAL_MATRIX_RESIZE(tmp2,3,3);
      MAL_S3x3_C_eq_A_by_B(tmp2_3d,aBody->R,aBody->getInertie()); 
      for(unsigned int i=0;i<3;++i)
	for(unsigned int j=0;j<3;++j)
	  tmp2(i,j) = tmp2_3d(i,j);

      MAL_C_eq_A_by_B(rightoperand,tmp2,pAngularJacobian);
      
      outjacobian = outjacobian + rightoperand;
    }
  
}
