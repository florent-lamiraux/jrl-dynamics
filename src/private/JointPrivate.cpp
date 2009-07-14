/*
  Copyright (c) 2005-2006, 
  @author Olivier Stasse, Oussama Kanoun, Fumio Kanehiro, Florent Lamiraux
   
  JRL-Japan, CNRS/AIST
 
  All rights reserved.

  Please refers to file License.txt for details on the license.
   
*/
#include "Debug.h"

#include "JointPrivate.h"
#include "DynamicBodyPrivate.h"


using namespace dynamicsJRLJapan;

JointPrivate::JointPrivate(int ltype, MAL_S3_VECTOR(,double) & laxe,
             float lquantite, MAL_S4x4_MATRIX(,double) & lpose):
  m_inGlobalFrame(false),
  m_type(ltype),
  m_axe(laxe),
  m_quantity(lquantite),
  m_poseInParentFrame(lpose),
  m_FatherJoint(0),
  m_IDinActuated(-1)
{
  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);
  CreateLimitsArray();
}

JointPrivate::JointPrivate(int ltype, MAL_S3_VECTOR(,double) & laxe,
             float lquantite, MAL_S3_VECTOR(,double) & translationStatic):
  m_inGlobalFrame(false),
  m_type(ltype),
  m_axe(laxe),
  m_quantity(lquantite),
  m_FatherJoint(0),
  m_IDinActuated(-1)
{
  MAL_S4x4_MATRIX_SET_IDENTITY(m_poseInParentFrame);
  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,0,3) = translationStatic[0];
  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,1,3) = translationStatic[1];
  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,2,3) = translationStatic[2];
  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);

  CreateLimitsArray();
}

JointPrivate::JointPrivate(int ltype, MAL_S3_VECTOR(,double) & laxe,
             float lquantite):
  m_inGlobalFrame(false),
  m_type(ltype),
  m_axe(laxe),
  m_quantity(lquantite),
  m_FatherJoint(0),
  m_IDinActuated(-1)
{
  MAL_S4x4_MATRIX_SET_IDENTITY(m_poseInParentFrame);
  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);
  CreateLimitsArray();
}

JointPrivate::JointPrivate(const JointPrivate &r)
{
  m_type = r.type();
  m_axe = r.axe();
  m_quantity=r.quantity();
  m_poseInParentFrame=r.pose();
  m_FatherJoint = 0;
  m_Name=r.getName();
  m_IDinActuated=r.getIDinActuated();
  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);
  m_inGlobalFrame=r.m_inGlobalFrame;
  CreateLimitsArray();

  for(unsigned int i=0;i<numberDof();i++)
    {
      m_LowerLimits[i] = r.lowerBound(i);
      m_UpperLimits[i] = r.upperBound(i);
      m_LowerVelocityLimits[i] = r.lowerVelocityBound(i);
      m_UpperVelocityLimits[i] = r.upperVelocityBound(i);
    }

}

JointPrivate::JointPrivate():
  m_inGlobalFrame(false),
  m_quantity(0.0),
  m_FatherJoint(0),
  m_IDinActuated(-1)
{
  MAL_S3_VECTOR_ACCESS(m_axe,0) = 0.0;
  MAL_S3_VECTOR_ACCESS(m_axe,1) = 0.0;
  MAL_S3_VECTOR_ACCESS(m_axe,2) = 0.0;
  MAL_S4x4_MATRIX_SET_IDENTITY(m_poseInParentFrame);

  m_type = FREE_JOINT;
  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);
  CreateLimitsArray();
}

JointPrivate::~JointPrivate()
{}

void JointPrivate::CreateLimitsArray()
{
  if (numberDof()!=0)
    {
      m_LowerLimits.resize(numberDof());
      m_UpperLimits.resize(numberDof());
      m_LowerVelocityLimits.resize(numberDof());
      m_UpperVelocityLimits.resize(numberDof());
      for (unsigned int i=0; i<numberDof(); i++)
        {
	  m_LowerLimits[i] = 0;
	  m_UpperLimits[i] = 0;
	  m_LowerVelocityLimits[i] = 0;
	  m_UpperVelocityLimits[i] = 0;
        }
    }
  else
    {
      m_LowerLimits.clear();
      m_UpperLimits.clear();
      m_LowerVelocityLimits.clear();
      m_UpperVelocityLimits.clear();
    }
}


void JointPrivate::computeLocalAndGlobalPose()
{
  if (m_FatherJoint==0)
    return;

  if (m_inGlobalFrame)
    {
      /*
	The pose of the joint has been defined in global frame at construction. 
	Compute pose in local frame of parent joint.
      */

      /* Get global pose of parent joint */
      MAL_S4x4_MATRIX(, double) invParentGlobalPose;
      MAL_S4x4_INVERSE(m_FatherJoint->m_globalPoseAtConstruction, invParentGlobalPose, double);
      MAL_S4x4_MATRIX(, double) jointGlobalPose = m_globalPoseAtConstruction;
      /*
	parent     /  global \  -1   global
	R         = | R        |     R
	joint      \  parent /       joint
      */

      m_poseInParentFrame = MAL_S4x4_RET_A_by_B(invParentGlobalPose, jointGlobalPose);
      ODEBUG2("m_FatherJoint->m_globalPoseAtConstruction=" << m_FatherJoint->m_globalPoseAtConstruction);
      ODEBUG2("invParentGlobalPose=" << invParentGlobalPose);
      ODEBUG2("jointGlobalPose=" << jointGlobalPose);
    }
  else
    {
      /*
	The pose of the joint has been defined in local frame of parent joint at construction.
	Compute pose in global frame.
      */

      /*
	global       global      global
	R         =  R           R
	joint        parent      joint
      */
      m_globalPoseAtConstruction = MAL_S4x4_RET_A_by_B(m_FatherJoint->m_globalPoseAtConstruction,
						       m_poseInParentFrame);

    }
}

JointPrivate & JointPrivate::operator=(const JointPrivate & r)
{
  m_type = r.type();
  m_axe = r.axe();
  m_quantity=r.quantity();
  m_poseInParentFrame=r.pose();
  m_Name = r.getName();
  m_IDinActuated = r.getIDinActuated();
  m_inGlobalFrame = r.m_inGlobalFrame;

  CreateLimitsArray();
  for(unsigned int i=0;i<numberDof();i++)
    {
      m_LowerLimits[i] = r.lowerBound(i);
      m_UpperLimits[i] = r.upperBound(i);
      m_LowerVelocityLimits[i] = r.lowerVelocityBound(i);
      m_UpperVelocityLimits[i] = r.upperVelocityBound(i);
    }
  return *this;
};


/***********************************************/
/* Implementation of the generic JRL interface */
/***********************************************/

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

std::vector<JointPrivate*> JointPrivate::jointsFromRootToThisJoint() const
{
  return m_FromRootToThisJoint;
}

const MAL_S4x4_MATRIX(,double) & JointPrivate::currentTransformation() const
{
  DynamicBodyPrivate *m_DBody = (DynamicBodyPrivate *) m_Body;
  return m_DBody->m_transformation;
}

CjrlRigidVelocity JointPrivate::jointVelocity()
{

  DynamicBodyPrivate *m_DBody = dynamic_cast<DynamicBodyPrivate *>(m_Body);
  CjrlRigidVelocity ajrlRV(m_DBody->v0,m_DBody->w);
  return ajrlRV;
}

void JointPrivate::computeSubTreeMCom()
{
  for (attId = 0; attId< 3;attId++)
    attSTmcom[attId] = linkedDBody()->massCoef()*linkedDBody()->w_c[attId];

  attSTcoef = linkedDBody()->massCoef();

  for (attId = 0; attId< countChildJoints();attId++)
    {
      m_Children[attId]->computeSubTreeMCom();
      attSTmcom += m_Children[attId]->subTreeMCom();
      attSTcoef += m_Children[attId]->subTreeCoef();
    }
}

void JointPrivate::computeSubTreeMComExceptChild(const CjrlJoint* inJoint)
{
  for (attId = 0; attId< 3;attId++)
    attSTmcom[attId] = linkedDBody()->massCoef()*linkedDBody()->w_c[attId];

  attSTcoef = linkedDBody()->massCoef();

  for (attId = 0; attId< countChildJoints();attId++)
    {
      if (inJoint == m_Children[attId])
	continue;
      m_Children[attId]->computeSubTreeMCom();
      attSTmcom += m_Children[attId]->subTreeMCom();
      attSTcoef += m_Children[attId]->subTreeCoef();
    }
}

void JointPrivate::subTreeMCom(const vector3d& inReplacement)
{
  attSTmcom = inReplacement;
}

const vector3d& JointPrivate::subTreeMCom() const
{
  return attSTmcom;
}

double JointPrivate::subTreeCoef()
{
  return attSTcoef;
}

void JointPrivate::subTreeCoef(double inReplacement)
{
  attSTcoef = inReplacement;
}

CjrlRigidAcceleration JointPrivate::jointAcceleration()
{
  // TODO : Update the member of this object
  // TODO : when calling ForwardDynamics.
  // TODO : This will avoid the dynamic cast.
  MAL_S3_VECTOR(,double) a,b;

  if (m_Body!=0)
    {
      DynamicBodyPrivate *m_DBody = dynamic_cast<DynamicBodyPrivate *>(m_Body);

      a = m_DBody->dv;
      b = m_DBody->dw;
    }
  CjrlRigidAcceleration ajrlRA(a,b);

  return ajrlRA;

}

unsigned int JointPrivate::numberDof() const
{
  unsigned int r=0;

  switch(m_type)
    {
    case (FREE_JOINT):
      r=6;
      break;
    case (FIX_JOINT):
      r=0;
      break;
    case (REVOLUTE_JOINT):
      r=1;
      break;
    case (PRISMATIC_JOINT):
      r=1;
      break;
    }
  return r;
}

const MAL_MATRIX(,double) & JointPrivate::jacobianJointWrtConfig() const
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
  
  ODEBUG("Size of the jacobian :" << m_FromRootToThis.size()-1);
  
  for(unsigned int i=0;i<m_FromRootToThisJoint.size();i++)
    {
      MAL_VECTOR_DIM(LinearAndAngularVelocity,double,6);

      DynamicBodyPrivate * aBody= m_FromRootToThisJoint[i]->linkedDBody();
      JointPrivate * aJoint = m_FromRootToThisJoint[i];
      
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
  if (outJ.size1() !=6 || outJ.size2() != m_J.size2())
    {
      outJ.resize(6,m_J.size2(),false);
    }
  outJ.clear();


  DynamicBodyPrivate * FinalBody = (DynamicBodyPrivate *)m_Body;

  vector3d pn = FinalBody->p + MAL_S3x3_RET_A_by_B(FinalBody->R, inPointJointFrame);
  getJacobianWorldPointWrtConfig(pn, outJ);
}


CjrlBody* JointPrivate::linkedBody() const
{
  return m_Body;
}

DynamicBodyPrivate* JointPrivate::linkedDBody() const
{
  return m_dynBody;
}

void JointPrivate::setLinkedBody(CjrlBody& inBody)
{
  m_Body = &inBody;
  m_dynBody = (DynamicBodyPrivate*)m_Body;
}

void JointPrivate::SetFatherJoint(JointPrivate *aFather)
{
  m_FatherJoint = aFather;

  m_FromRootToThis.clear();
  m_FromRootToThisJoint.clear();

  m_FromRootToThis.push_back(this);
  m_FromRootToThisJoint.push_back(this);

  CjrlJoint* aJoint = m_FatherJoint;
  while(aJoint!=0)
    {
      m_FromRootToThis.insert(m_FromRootToThis.begin(),aJoint);
      m_FromRootToThisJoint.insert(m_FromRootToThisJoint.begin(),(JointPrivate*)aJoint);
      aJoint = aJoint->parentJoint();
    }
  computeLocalAndGlobalPose();
}

const MAL_S4x4_MATRIX(,double) & JointPrivate::initialPosition()
{
  if (m_Body!=0)
    {
      DynamicBodyPrivate *aDB = (DynamicBodyPrivate *) m_Body;
      ODEBUG("JointPrivate Name " << m_Name << " " << m_Body);
      for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,i,j) = aDB->R(i,j);
      for(int i=0;i<3;i++)
	MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,i,3) = aDB->p(i);
      ODEBUG( m_poseInParentFrame);

    }
  return m_poseInParentFrame;
}

void JointPrivate::UpdatePoseFrom6DOFsVector(MAL_VECTOR(,double) a6DVector)
{
  // Update the orientation of the joint.
  // Takes the three euler joints

  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,0,3) = a6DVector(0);
  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,1,3) = a6DVector(1);
  MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,2,3) = a6DVector(2);

  DynamicBodyPrivate* body = dynamic_cast<DynamicBodyPrivate*>(m_Body);
  if (!body)
    {
      std::cerr << "m_Body is not an instance of DynamicBodyPrivate" << std::endl;
    }
  body->p[0] = a6DVector(0);
  body->p[1] = a6DVector(1);
  body->p[2] = a6DVector(2);

  MAL_S3x3_MATRIX(,double) D,B,C,A;
  double CosTheta, SinTheta,
    CosPhi, SinPhi,
    CosPsi, SinPsi;


  CosPsi = cos(a6DVector(3));
  SinPsi = sin(a6DVector(3));
  CosTheta = cos(a6DVector(4));
  SinTheta = sin(a6DVector(4));
  CosPhi = cos(a6DVector(5));
  SinPhi = sin(a6DVector(5));

  /*
    D(0,0) =       1; D(0,1) =        0; D(0,2) = 0;
    D(1,0) =       0; D(1,1) =   CosPsi; D(1,2) = -SinPsi;
    D(2,0) =       0; D(2,1) =   SinPsi; D(2,2) = CosPsi;

    C(0,0) =  CosTheta; C(0,1) =        0; C(0,2) = SinTheta;
    C(1,0) =         0; C(1,1) =        1; C(1,2) = 0;
    C(2,0) = -SinTheta; C(2,1) =        0; C(2,2) = CosTheta;

    B(0,0) =  CosPhi; B(0,1) = -SinPhi; B(0,2) = 0;
    B(1,0) =  SinPhi; B(1,1) = CosPhi;  B(1,2) = 0;
    B(2,0) =       0; B(2,1) =      0;  B(2,2) = 1;

    MAL_S3x3_MATRIX(,double) tmp;
    MAL_S3x3_C_eq_A_by_B(tmp,C,D);
    MAL_S3x3_C_eq_A_by_B(A,B,tmp);

    body->R = A;
  */

  //Formulae for the above commented rotation composition
  A(0,0) = CosTheta * CosPhi ;
  A(1,0) = CosTheta * SinPhi;
  A(2,0) = -SinTheta;

  A(0,1) = CosPhi * SinPsi * SinTheta - CosPsi * SinPhi;
  A(1,1) = CosPsi * CosPhi + SinPsi * SinTheta * SinPhi;
  A(2,1) = CosTheta * SinPsi;

  A(0,2) = CosPsi * CosPhi * SinTheta + SinPhi * SinPsi;
  A(1,2) = - CosPhi * SinPsi + CosPsi * SinTheta * SinPhi;
  A(2,2) = CosPsi * CosTheta;

  body->R = A;

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      MAL_S4x4_MATRIX_ACCESS_I_J(m_poseInParentFrame,i,j) = A(i,j);

  ODEBUG("m_poseInParentFrame : " << m_poseInParentFrame <<
	 " A: "<<endl << A <<
	 " tmp " << endl << tmp <<
	 "C " << endl << C <<
	 "D " << endl << D <<
	 "B " << endl << B );
  body->m_transformation = m_poseInParentFrame;
}

void JointPrivate::UpdateVelocityFrom2x3DOFsVector(MAL_S3_VECTOR(,double) & aLinearVelocity,
					    MAL_S3_VECTOR(,double) & anAngularVelocity)
{
  m_RigidVelocity.linearVelocity(aLinearVelocity);
  m_RigidVelocity.rotationVelocity(anAngularVelocity);
}


void JointPrivate::RodriguesRotation(vector3d& inAxis, double inAngle, matrix3d& outRotation)
{
  double norm_w = MAL_S3_VECTOR_NORM(inAxis);
  if (norm_w < 10e-7)
    {
      MAL_S3x3_MATRIX_SET_IDENTITY(outRotation);
    }
  else
    {
      double th = norm_w * inAngle;
      wn3d = inAxis / norm_w;
      double ct = cos(th);
      double lct= (1-ct);
      double st = sin(th);
      outRotation(0,0) = ct + wn3d[0]*wn3d[0]* lct;
      outRotation(0,1) = wn3d[0]*wn3d[1]*lct-wn3d[2]*st;
      outRotation(0,2) = wn3d[1] * st+wn3d[0]*wn3d[2]*lct;
      outRotation(1,0) = wn3d[2]*st +wn3d[0]*wn3d[1]*lct;
      outRotation(1,1) = ct + wn3d[1]*wn3d[1]*lct;
      outRotation(1,2) = -wn3d[0]*st+wn3d[1]*wn3d[2]*lct;
      outRotation(2,0) = -wn3d[1]*st+wn3d[0]*wn3d[2]*lct;
      outRotation(2,1) = wn3d[0]*st + wn3d[1]*wn3d[2]*lct;
      outRotation(2,2) = ct + wn3d[2]*wn3d[2]*lct;
    }
}

JointFreeflyerPrivate::JointFreeflyerPrivate(const MAL_S4x4_MATRIX(,double) &inInitialPosition)
{
  type(JointPrivate::FREE_JOINT);
  m_inGlobalFrame = true;
  m_globalPoseAtConstruction = inInitialPosition;

  ODEBUG2("freeflyer: inInitialPosition" << inInitialPosition);
  MAL_S3_VECTOR(axis, double);

  MAL_S3_VECTOR_ACCESS(axis,0) = 1.0;
  MAL_S3_VECTOR_ACCESS(axis,1) = 0.0;
  MAL_S3_VECTOR_ACCESS(axis,2) = 0.0;

  axe(axis);
}

JointFreeflyerPrivate::~JointFreeflyerPrivate()
{
}

JointAnchorPrivate::JointAnchorPrivate(const MAL_S4x4_MATRIX(,double) &inInitialPosition)
{
  type(JointPrivate::FIX_JOINT);
  m_inGlobalFrame = true;
  m_globalPoseAtConstruction = inInitialPosition;

  ODEBUG2("anchor: inInitialPosition" << inInitialPosition);
  MAL_S3_VECTOR(axis, double);

  MAL_S3_VECTOR_ACCESS(axis,0) = 1.0;
  MAL_S3_VECTOR_ACCESS(axis,1) = 0.0;
  MAL_S3_VECTOR_ACCESS(axis,2) = 0.0;

  axe(axis);
}

JointAnchorPrivate::~JointAnchorPrivate()
{
}

JointRotationPrivate::JointRotationPrivate(const MAL_S4x4_MATRIX(,double) &inInitialPosition)
{
  type(JointPrivate::REVOLUTE_JOINT);
  m_inGlobalFrame = true;
  m_globalPoseAtConstruction = inInitialPosition;

  ODEBUG2("rotation: inInitialPosition" << inInitialPosition);

  MAL_S3_VECTOR(axis, double);

  MAL_S3_VECTOR_ACCESS(axis,0) = 1.0;
  MAL_S3_VECTOR_ACCESS(axis,1) = 0.0;
  MAL_S3_VECTOR_ACCESS(axis,2) = 0.0;

  axe(axis);
}

JointRotationPrivate::~JointRotationPrivate()
{
}

JointTranslationPrivate::JointTranslationPrivate(const MAL_S4x4_MATRIX(,double) &inInitialPosition)
{
  type(JointPrivate::PRISMATIC_JOINT);
  m_inGlobalFrame = true;
  m_globalPoseAtConstruction = inInitialPosition;

  ODEBUG2("translation: inInitialPosition" << inInitialPosition);
  MAL_S3_VECTOR(axis, double);

  MAL_S3_VECTOR_ACCESS(axis,0) = 1.0;
  MAL_S3_VECTOR_ACCESS(axis,1) = 0.0;
  MAL_S3_VECTOR_ACCESS(axis,2) = 0.0;

  axe(axis);
}

JointTranslationPrivate::~JointTranslationPrivate()
{}


bool JointPrivate::updateTransformation(const vectorN& inDofVector)
{
  if (rankInConfiguration() > inDofVector.size() -1 )
    {
      std::cout << "JointTranslationPrivate::updateTransformation(). Inappropriate configuration vector.\n";
      return false;
    }

  switch (type())
    {
    case JointPrivate::REVOLUTE_JOINT :
      {
	DynamicBodyPrivate* body = (DynamicBodyPrivate*)(linkedBody());
	DynamicBodyPrivate* parentbody = (DynamicBodyPrivate*)(parentJoint()->linkedBody());

	body->q = inDofVector(rankInConfiguration());
	quantity( body->q);

	RodriguesRotation(body->a, body->q, localR);

	MAL_S3x3_MATRIX(,double) Rtmp;
	MAL_S3x3_C_eq_A_by_B(Rtmp ,parentbody->R , body->R_static);
	MAL_S3x3_C_eq_A_by_B(body->R , Rtmp, localR);
	body->p = parentbody->p + MAL_S3x3_RET_A_by_B(parentbody->R,body->b);

	for( unsigned int i=0;i<3;i++)
	  for(unsigned int j=0;j<3;j++)
	    MAL_S4x4_MATRIX_ACCESS_I_J(body->m_transformation,i,j) = body->R(i,j);

	for( unsigned int i=0;i<3;i++)
	  MAL_S4x4_MATRIX_ACCESS_I_J(body->m_transformation,i,3) = body->p(i);
      }
      break;

    case JointPrivate::FREE_JOINT :
      {
	if (dof6D.size() != 6)
	  dof6D.resize(6,false);

	for (unsigned int i=0; i<6; i++)
	  dof6D(i) = inDofVector(rankInConfiguration() + i);

	UpdatePoseFrom6DOFsVector(dof6D);
      }
      break;

    case JointPrivate::PRISMATIC_JOINT :
      {
	DynamicBodyPrivate* body = (DynamicBodyPrivate*)(linkedBody());
	DynamicBodyPrivate* parentbody = (DynamicBodyPrivate*)(parentJoint()->linkedBody());

	body->q = inDofVector(rankInConfiguration());
	quantity( body->q);

	for (unsigned int i = 0; i<3; i++)
	  vek[i] *= body->q;

	body->R = parentbody->R;

	MAL_S3x3_C_eq_A_by_B(wn3d, body->R, vek);

	body->p = parentbody->p+ MAL_S3x3_RET_A_by_B(parentbody->R,body->b) + wn3d;
      }
      break;
    }
  return true;
}