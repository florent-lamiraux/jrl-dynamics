/*
 * Copyright 2010,
 *
 * Oussama Kanoun
 * Florent Lamiraux
 * Olivier Stasse
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
/*! Computation of the dynamic aspect for a humanoid robot.
*/

#ifndef HUMDYNMULTIBODYPRIVATE_H
#define HUMDYNMULTIBODYPRIVATE_H

#include <vector>

#include "jrl/dynamics/foot.hh"

#include "jrl/mal/matrixabstractlayer.hh"
#include "DynMultiBodyPrivate.h"

#include <abstract-robot-dynamics/humanoid-dynamic-robot.hh>
namespace dynamicsJRLJapan
{
  class HumanoidSpecificities;

  struct HumanoidDMBComparison
  {
    bool operator()(int a,int b) const
    {
      return a < b;
    }
  };

  /** \brief This class implements the functionnalities specific to a dynamic model for a humanoid.

  This includes a direct access to the joints reprensenting the hands, the foot, and the gaze.
  This specific class is the specialization of the generic class CjrlHumanoidDynamicRobot.

  */
  class HumDynMultiBodyPrivate:     public virtual CjrlHumanoidDynamicRobot,
  public DynMultiBodyPrivate
    {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    private:

      /** \brief Store the Left Wrist Joint */
      CjrlJoint * m_LeftWristJoint;

      /** \brief Store the Right Wrist Joint */
      CjrlJoint * m_RightWristJoint;

      /** \brief Store the Left Foot */
      CjrlFoot * m_LeftFoot;

      /** \brief Store the Right Foot */
      CjrlFoot *m_RightFoot;

      /** \brief Store the Left Ankle Joint */
      CjrlJoint* m_LeftAnkleJoint;

      /** \brief Store the Right Ankle Joint */
      CjrlJoint* m_RightAnkleJoint;

      /** \brief Store the Left hand */
      CjrlHand * m_leftHand;

      /** \brief Store the Right hand */
      CjrlHand * m_rightHand;

      /** \brief Store the Gaze Joint */
      CjrlJoint * m_GazeJoint;

      /** \brief Store the Waist Joint */
      CjrlJoint * m_WaistJoint;

      /*! \brief Store the chest. */
      CjrlJoint * m_ChestJoint;

      /** \name Gaze related store */

      /** \brief Set the direction of the line. */
      vector3d m_LineVector;

      /** \brief Set the point through which the line is going. */
      vector3d m_LinePoint;

      /** @} */

      /*! Object to store the specificities to an instance of a humanoid robot. */
      HumanoidSpecificities *m_HS;

      /*! Distance between the ankle and the soil. */
      double m_AnkleSoilDistance;

      /*! Distance between some axis of the hip. */
      vector3d m_Dt;

      /*! \name Parameters related to the distance
	between the CoM and the hip.
	@{
      */

      /*! Static translation from the CoM to the left
       hip. */
      vector3d m_StaticToTheLeftHip,
	m_StaticToTheRightHip;

      /*! Static translation from the CoM to the right
       hip. */
      vector3d m_TranslationToTheRightHip,
	m_TranslationToTheLeftHip;

      vector3d m_ZeroMomentumPoint;
      /*! @} */

    public:

      /** \name Constructors and destructors
	  @{
      */
      /*! Constructor */

      /*! Default constructor: assume that a file name HumanoidSpecificities.xml
	 is present in the current directory, and creates a Dynamic Multi Body
	 object instead of a more generic multibody.
	 Good for use by default.
      */
      HumDynMultiBodyPrivate(void);

      /*! Advanced constructor: */
      HumDynMultiBodyPrivate(const DynMultiBodyPrivate& inDynamicMultiBody,
			       string aFileNameForHumanoidSpecificities);

      /*! Destructor */
      virtual ~HumDynMultiBodyPrivate();

      /** @} */

      /*! This method creates the link between the Humanoid Specificities
	 object and the Dynamic MultiBody fields. */
      void LinkBetweenJointsAndEndEffectorSemantic();

      /*! Set Humanoid Specificties file. */
      void SetHumanoidSpecificitiesFile(string &aFileName);

      /*! Get pointer on the information specific to the humanoid */
      inline HumanoidSpecificities * getHumanoidSpecificities() const
	{return m_HS;};


      /** \name jrlHumanoidDynamicRobot Interface */

      /**
	 \name Joints specific to humanoid robots
      */

      /**
	\brief Set the pointer to the Waist joint
      */
      void waist(CjrlJoint* inWaist);

      /**
	 \brief Get a pointer to the Waist
      */
      CjrlJoint* waist() const;

      /**
      \brief Set the pointer to the left Wrist joint.
      */
      void leftWrist(CjrlJoint *inLeftWrist);


      /**
      \brief Get a pointer to the left Wrist.
      */
      CjrlJoint *leftWrist() const;

      /**
      \brief Set the pointer to the right Wrist joint.
      */
      void rightWrist(CjrlJoint *inRightWrist);

      /**
      \brief Get a pointer to the right Wrist.
      */
      CjrlJoint *rightWrist() const;


      /**
        \brief Set the pointer to the right hand
      */
      virtual void rightHand(CjrlHand* inRightHand);

      /**
        \brief Get a pointer to the right hand
      */
      virtual CjrlHand* rightHand() const;

      /**
        \brief Set the pointer to the left hand
      */
      virtual void leftHand(CjrlHand* inLeftHand);

      /**
        \brief Get a pointer to the left hand
      */
      virtual CjrlHand* leftHand() const;

      /**
	 \brief Get the hand clench value. This is a scalar value ranging between 0
	 and 1 which describes the hand clench (0 for open and 1 for closed hand)
      */
      virtual double getHandClench(CjrlHand* inHand);

      /**
	 \brief Set the hand clench value. This is a scalar value ranging
	 between 0 and 1 which describes the hand clench (0 for open and 1 for closed hand).
      */
      virtual bool  setHandClench(CjrlHand* inHand, double inClenchingValue);

      /**
	 \brief Set the pointer to the left ankle joint.
      */
      void leftAnkle(CjrlJoint* inLefAnkle);

      /**
	 \brief Get a pointer to the left ankle.
      */
      CjrlJoint* leftAnkle() const;

      /**
	 \brief Set the pointer to the right ankle joint.
      */
      void rightAnkle(CjrlJoint* inRightAnkle);

      /**
	 \brief Get a pointer to the right ankle.
      */
      CjrlJoint* rightAnkle() const;

      /**
	 \brief Set the pointer to the left foot joint.
      */
      void leftFoot(CjrlFoot *inLeftFoot);

      /**
	  \brief Get a pointer to the left foot.
      */
      CjrlFoot *leftFoot() const;

      /**
	 \brief Set the pointer to the right foot joint.
      */
      void rightFoot(CjrlFoot *inRightFoot);

      /**
	  \brief Get a pointer to the right foot.
      */
      CjrlFoot *rightFoot() const;

      /**
	  \brief Set gaze joint

	  \note  For most humanoid robots, the gaze joint is the head.
      */
      inline void gazeJoint(CjrlJoint *inGazeJoint)
	{ m_GazeJoint = (JointPrivate *)inGazeJoint; }

      /**
	 \brief Get gaze joint
      */
      CjrlJoint *gazeJoint() const
	{ return m_GazeJoint; }

      /**
	 \brief Set the gaze in the local frame of the gaze joint.

	 \note The gaze is defined as a straight line linked to the gaze joint.
	 @param inVector: A 3D vector which define the direction of the gaze.
	 @param inPoint: A 3D point by which the line of direction \a inVector
	 goes through.

	 Those two paramaters defines a line in the head reference frame defining
	 a gaze direction.

      */
      inline void gaze( const vector3d& inVector, const vector3d & inPoint)
	{ m_LineVector = inVector; m_LinePoint = inPoint;};

      /**
      \brief Get a point on the gaze straight line
      */
      const vector3d & gazeOrigin() const {return m_LinePoint;}

      /**
	 \brief Get the direction of gaze
      */
      const vector3d & gazeDirection() const {return m_LineVector;}


      /*! \brief Specify which joint is the chest. */
      void chest(CjrlJoint *);

      /*! \brief Returns joint which is the first joint of the chest. */
      CjrlJoint * chest() const;

      /**
	 \@}
      */

      /** \name Methods related to fixed joints.
	  @{
      */

      /**
	 \@}
      */

      /**
	 \name Zero momentum point
      */

      /**
	 \brief Compute the coordinates of the Zero Momentum Point.
      */
      const vector3d & zeroMomentumPoint() const;

      /**
	  \brief Trigger the computation of the Zero Momentum Point.
      */
      void ComputingZeroMomentumPoint();

      /**
	 @}
      */

      /**
      @}
      */

      /**
	  \name Forward kinematics and dynamics
      */


      /**
	 \brief Compute forward kinematics.

	 Update the position, velocity and accelerations of each
	 joint wrt \f${\bf {q}}\f$, \f${\bf \dot{q}}\f$, \f${\bf \ddot{q}}\f$.

	 Computes the ZMP position.

      */
      bool computeForwardKinematics() ;

      /**
	 @}
      */

      /**
	  \name Jacobian fonctions
      */

      /**
	 \brief Get the jacobian of a joint wrt to internal configuration variables assuming a joint is fixed.

	 Fixed joint is first fixed joint in vector.
	 \return true if there is at least one fixed joint, false otherwise.
      */
      bool jacobianJointWrtFixedJoint(CjrlJoint *inJoint,
				      matrixNxP & outJacobian);

      /**
	 @}
      */

    };

}
#endif /* HUMDYNMULTIBODYPRIVATE_H */
