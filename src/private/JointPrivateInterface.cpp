/*
 * Copyright 2010,
 *
 * Olivier Stasse
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

/* \file This part implements the generic JRL interface
*/

#include "Debug.h"

#include <jrl/mal/matrixabstractlayer.hh>
#include <jrl/dynamics/dynamicbody.hh>
#include "JointPrivate.h"
#include "DynamicBodyPrivate.h"

using namespace dynamicsJRLJapan;

CjrlJoint* JointPrivate::parentJoint() const
{
  return m_FatherJoint;
}

bool JointPrivate::addChildJoint(CjrlJoint& aJoint)
{
  JointPrivate * pjoint = (JointPrivate *)&aJoint;
  for(unsigned int li =0; li < m_Children.size();li++)
    {
      if (m_Children[li]==pjoint)
	return true;
    }
  // Make sure I went through this part.
  ODEBUG("Set father joint : " << pjoint->getName()
	 << " "  << getName());
  pjoint->SetFatherJoint(this);
  m_Children.push_back(pjoint);
  return true;
}

unsigned int JointPrivate::countChildJoints() const
{
  return m_Children.size();
}

CjrlJoint* JointPrivate::childJoint(unsigned int givenRank) const
{
  if (givenRank<m_Children.size())
    return m_Children[givenRank];

  return 0;
}

JointPrivate* JointPrivate::child_JointPrivate(unsigned int givenRank) const
{
  if (givenRank<m_Children.size())
    return m_Children[givenRank];

  return 0;
}

std::vector<CjrlJoint*> JointPrivate::jointsFromRootToThis() const
{
  return m_FromRootToThis;
}


const MAL_S4x4_MATRIX_TYPE(double) & JointPrivate::currentTransformation() const
{
  if (m_dynBody==0)
    return m_globalPoseAtConstruction;
  return m_dynBody->m_transformation;
}

CjrlRigidVelocity JointPrivate::jointVelocity() const
{

  CjrlRigidVelocity ajrlRV(m_dynBody->v0,m_dynBody->w);
  return ajrlRV;
}

void JointPrivate::computeSubTreeMCom()
{

  for (unsigned int Id = 0; Id< 3;Id++)
    m_STmcom[Id] = m_dynBody->massCoef()*m_dynBody->w_c[Id];

  m_STcoef = m_dynBody->massCoef();

  for (unsigned int Id = 0; Id< countChildJoints();Id++)
    {
      m_Children[Id]->computeSubTreeMCom();
      m_STmcom += m_Children[Id]->subTreeMCom();
      m_STcoef += m_Children[Id]->subTreeCoef();
    }
}

void JointPrivate::computeSubTreeMComExceptChild(const CjrlJoint* inJoint)
{
  for (unsigned int Id = 0; Id< 3;Id++)
    m_STmcom[Id] = m_dynBody->massCoef()*m_dynBody->w_c[Id];

  m_STcoef = m_dynBody->massCoef();

  for (unsigned int Id = 0; Id< countChildJoints();Id++)
    {
      if (inJoint == m_Children[Id])
	continue;
      m_Children[Id]->computeSubTreeMCom();
      m_STmcom += m_Children[Id]->subTreeMCom();
      m_STcoef += m_Children[Id]->subTreeCoef();
    }
}

void JointPrivate::subTreeMCom(const vector3d& inReplacement)
{
  m_STmcom = inReplacement;
}

const vector3d& JointPrivate::subTreeMCom() const
{
  return m_STmcom;
}

double JointPrivate::subTreeCoef()
{
  return m_STcoef;
}

void JointPrivate::subTreeCoef(double inReplacement)
{
  m_STcoef = inReplacement;
}

CjrlRigidAcceleration JointPrivate::jointAcceleration() const
{
  MAL_S3_VECTOR_TYPE(double) a,b;

  if (m_dynBody!=0)
    {
      a = m_dynBody->dv;
      b = m_dynBody->dw;
    }
  CjrlRigidAcceleration ajrlRA(a,b);

  return ajrlRA;

}


const MAL_MATRIX_TYPE(double) & JointPrivate::jacobianJointWrtConfig() const
{
  return m_J;
}

void JointPrivate::resizeJacobianJointWrtConfig(int lNbDofs)
{
  MAL_MATRIX_RESIZE(m_J,6,lNbDofs);
  MAL_MATRIX_FILL(m_J,0.0);

}


void JointPrivate::computeJacobianJointWrtConfig()
{
  DynamicBodyPrivate * FinalBody = (DynamicBodyPrivate *)m_Body;
  getJacobianWorldPointWrtConfig(FinalBody->p, m_J);
}

void JointPrivate::getJacobianWorldPointWrtConfig(const vector3d& inPointWorldFrame,
						  matrixNxP& outJ) const
{
  vector3d dp,lv;

  ODEBUG("Size of the jacobian :" << m_FromRootToThisPrivate.size()-1);

  for(unsigned int i=0;i<m_FromRootToThisPrivate.size();i++)
    {
      MAL_VECTOR_DIM(LinearAndAngularVelocity,double,6);

      DynamicBodyPrivate * aBody= m_FromRootToThisPrivate[i]->m_dynBody;
      JointPrivate * aJoint = m_FromRootToThisPrivate[i];

      unsigned int lcol = aJoint->stateVectorPosition();
      ODEBUG("JointPrivate: " << aJoint->getName() << " " << lcol);
      dp = inPointWorldFrame - aBody->p;

      MAL_S3_VECTOR_CROSS_PRODUCT(lv,aBody->w_a,dp);

      switch (aJoint->type())
        {

        case JointPrivate::REVOLUTE_JOINT:
	  for(int j=0;j<3;j++)
            {
	      outJ(j,lcol) =  lv[j];
	      outJ(j+3,lcol) = aBody->w_a[j];
            }
	  break;
        case JointPrivate::PRISMATIC_JOINT:
	  for(int j=0;j<3;j++)
            {
	      outJ(j,lcol) =  aBody->w_a[j];
	      outJ(j+3,lcol) = 0;
            }
	  break;
        case JointPrivate::FREE_JOINT:
	  //J =  I M = J11 J12
	  //     0 I   J21 J22
	  //
	  // with M = d(w x dp)/dw
	  //
	  for(int j=0;j<3;j++)
            {
	      for(int k=0;k<3;k++)
                {
		  // Computation of J11, J12 and J21
		  if (j!=k)
                    {
		      outJ(     j, lcol + k) =0.0;
		      outJ( j + 3, lcol + k + 3) = 0.0;
                    }
		  else
                    {
		      outJ(     j, lcol + k ) =1.0;
		      outJ( j + 3, lcol + k + 3) = 1.0;
                    }
		  outJ(j+3,k) = 0.0;
                }
            }
	  // Compute M
	  outJ( 0, lcol + 3 ) =      0;
	  outJ( 0 , lcol + 4 ) =   dp(2);
	  outJ( 0 , lcol + 5 ) = -dp(1);
	  outJ( 1, lcol + 3 ) = -dp(2);
	  outJ( 1 , lcol + 4 ) =      0 ;
	  outJ( 1 , lcol + 5 ) =  dp(0);
	  outJ( 2, lcol + 3 ) =  dp(1);
	  outJ( 2 , lcol + 4 ) =  -dp(0);
	  outJ( 2 , lcol + 5 ) =      0;
	  break;
        }
    }
}

/**
   \brief Get the jacobian of the point specified in local frame by inPointJointFrame.
   The output matrix outjacobian is automatically resized if necessary

*/
void JointPrivate::getJacobianPointWrtConfig(const vector3d& inPointJointFrame, matrixNxP& outJ) const
{
  if (MAL_MATRIX_NB_ROWS(outJ) !=6 || MAL_MATRIX_NB_COLS(outJ) != MAL_MATRIX_NB_COLS(m_J))
    {
      MAL_MATRIX_RESIZE(outJ,6,MAL_MATRIX_NB_COLS(m_J));
    }
  MAL_MATRIX_CLEAR(outJ);


  DynamicBodyPrivate * FinalBody = (DynamicBodyPrivate *)m_Body;

  vector3d pn = FinalBody->p + MAL_S3x3_RET_A_by_B(FinalBody->R, inPointJointFrame);
  getJacobianWorldPointWrtConfig(pn, outJ);
}


CjrlBody* JointPrivate::linkedBody() const
{
  return m_dynBody;
}

DynamicBodyPrivate* JointPrivate::linkedDBody() const
{
  return m_dynBody;
}

void JointPrivate::setLinkedBody(CjrlBody& inBody)
{
  DynamicBody* dynamicBody = dynamic_cast <DynamicBody*> (&inBody);
  if (dynamicBody) {
    m_Body = dynamicBody->m_privateObj;
    m_dynBody = dynamicBody->m_privateObj;
  } else {
    DynamicBodyPrivate* priv = dynamic_cast <DynamicBodyPrivate*> (&inBody);
    if (priv) {
      m_Body = priv;
      m_dynBody = priv;
    } else {
      abort ();
    }
  }
  resizeSpatialFields();
}

void JointPrivate::setLinkedDBody(DynamicBodyPrivate * inBody)
{
  m_dynBody = inBody;
  m_Body = m_dynBody;
  resizeSpatialFields();
}

void JointPrivate::SetFatherJoint(JointPrivate *aFather)
{
  m_FatherJoint = aFather;

  m_FromRootToThis.clear();
  m_FromRootToThisPrivate.clear();

  m_FromRootToThis.push_back(this);
  m_FromRootToThisPrivate.push_back(this);

  JointPrivate* aJoint = m_FatherJoint;
  while(aJoint!=0)
    {
      m_FromRootToThis.insert(m_FromRootToThis.begin(),aJoint);
      m_FromRootToThisPrivate.insert(m_FromRootToThisPrivate.begin(),aJoint);
      aJoint = aJoint->m_FatherJoint;
    }
}

const MAL_S4x4_MATRIX_TYPE(double) & JointPrivate::initialPosition() const
{
  return m_globalPoseAtConstructionNormalized;
}

