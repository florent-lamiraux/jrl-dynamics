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

void DynMultiBodyPrivate::BackwardDynamics(DynamicBodyPrivate & CurrentBody )
{
  ODEBUG("=====================");
  ODEBUG("Body : " << CurrentBody.getName());
  MAL_S3x3_MATRIX(,double) aRt;

  MAL_S3x3_MATRIX(,double) currentBodyR;
  currentBodyR = MAL_S3x3_RET_TRANSPOSE(CurrentBody.R);

  MAL_S3_VECTOR(,double) lg;
  lg(0) = 0.0;
  lg(1) = 0.0;
  lg(2) = -9.81;

  /* Compute the torque
   * with eq. (7.147) Spong RMC p. 277
   *
   *
   */
  MAL_S3_VECTOR(,double) firstterm,
    sndterm, thirdterm, fifthterm,tmp,tmp2;
  // There is no fourth term because it is the angular acceleration.

  /* Constant part */
  tmp = CurrentBody.dv_c - lg;
  CurrentBody.m_Force =  tmp * CurrentBody.mass();
  ODEBUG("tmp: " << tmp << " tmp2: " << tmp2 );
  /* 2nd term : -f_i x r_{i,ci} */
  vector3d lc = CurrentBody.localCenterOfMass();
  MAL_S3_VECTOR_CROSS_PRODUCT(sndterm,CurrentBody.m_Force, lc);


  /* 5th term : w_i x (I_i w_i)*/
  MAL_S3x3_MATRIX(,double) lI = CurrentBody.getInertie();
  tmp = MAL_S3x3_RET_A_by_B(lI,CurrentBody.w);
  MAL_S3_VECTOR_CROSS_PRODUCT(fifthterm,CurrentBody.w,tmp);

  CurrentBody.m_Torque =  CurrentBody.dw + fifthterm -sndterm;

  /* Compute with the force
   * eq. (7.146) Spong RMC p. 277
   * fi = R^i_{i+1} * f_{i+1} + m_i * a_{c,i} - m_i * g_i
   * g_i is the gravity express in the i reference frame.
   */


  int IndexChild = CurrentBody.child;
  ODEBUG( " Force from Acceleration + gravity: " << CurrentBody.m_Force <<
	   " Acceleration: " << CurrentBody.dv_c );
  while(IndexChild!=-1)
    {
      DynamicBodyPrivate *Child = m_listOfBodies[IndexChild];
      ODEBUG( "Child Bodies : " << Child->getName() );
      aRt = Child->Riip1;
      ODEBUG( "Riip1: " << aRt );
      // /* Force computation. */
      //// Other immediate child are sisters of the other immediate childs.
      ODEBUG("Child Force: " << Child->m_Force );
      tmp= MAL_S3x3_RET_A_by_B(aRt, Child->m_Force);
      ODEBUG("Riip1 * fi+1: " << tmp );
      CurrentBody.m_Force += tmp;

      /* Torque computation. */
      /* 1st term : R^i_{i+1} t_{i+1} */
      firstterm = MAL_S3x3_RET_A_by_B(aRt, Child->m_Torque);

      /* 3rd term : R_i_{i+1} f_{i+1} */
      vector3d rip1ci = Child->p - CurrentBody.w_c;
      ODEBUG( "rip1ci:" << rip1ci);
      MAL_S3_VECTOR_CROSS_PRODUCT(thirdterm,tmp, rip1ci);

      CurrentBody.m_Torque += firstterm + thirdterm;

      /* Body selection. */
      IndexChild = m_listOfBodies[IndexChild]->sister;
      if (IndexChild!=-1)
	Child=m_listOfBodies[IndexChild];
    }
  ODEBUG("Force for Body " << CurrentBody.getName() << " " << CurrentBody.m_Force );
  ODEBUG("=====================");
  // Update the vector related to the computed quantities.
  for(unsigned int i=0;i<m_StateVectorToJoint.size();i++)
    {
      unsigned int StateRankComputed=false;

      JointPrivate * aJoint = (JointPrivate *)m_JointVector[m_StateVectorToJoint[i]];
      if (aJoint!=0)
        {
	  DynamicBodyPrivate *aDB = (DynamicBodyPrivate *) aJoint->linkedBody();
	  if (aDB!=0)
            {
	      StateRankComputed = true;
	      for(unsigned int k=0;k<3;k++)
                {
		  m_Forces(i,k)=aDB->m_Force[k];
		  m_Torques(i,k) = aDB->m_Torque[k];
                }

            }
        }

      if (!StateRankComputed)
        {
	  for(unsigned int k=0;k<3;k++)
            {
	      m_Forces(i,k)=0.0;
	      m_Torques(i,k)=0.0;
            }
        }
    }
}
