/* @doc Object used to handle a foot 

   Copyright (c) 2009, 

   @author : 
   Olivier Stasse.

   JRL-Japan, CNRS/AIST

   All rights reserved.
   
   Please refers to file License.txt for details on the license.

*/

#ifndef _DYN_JRL_JAPAN_FOOT_H_
#define _DYN_JRL_JAPAN_FOOT_H_
#include <MatrixAbstractLayer/MatrixAbstractLayer.h>

#include "robotDynamics/jrlFoot.h"
#include "dynamicsJRLJapan/dll.h"

namespace dynamicsJRLJapan
{
  /*! \ingroup userclasses 
    This class represents a foot of a humanoid robot.
    It assumes some geometrical information available.
    They are described in more details in the original 
    class jrlFoot. */
  class DYN_JRL_JAPAN_EXPORT Foot: public CjrlFoot
  {
  public:

    /*! \brief Default constructor */
    Foot();
    Foot(const Foot &inFoot);

        /**
    \brief Destructor
     */
    virtual ~Foot();

    /*! Returns associated ankle. */
    virtual const CjrlJoint * associatedAnkle() const;

    /*! Returns associated ankle. */
    void setAssociatedAnkle(const CjrlJoint * inAssociatedAnkle);
    
    /** 
	\brief Get size of the rectagular sole
	
	\retval outLength length of the sole (see Figure)
	\retval outWidth width of the sole (see Figure)
	
    */
    virtual void getSoleSize(double &outLength, double &outWidth) const;

    /** 
	\brief Set size of the rectagular sole

	\param inLength length of the sole (see Figure)
	\param inWidth width of the sole (see Figure)

    */
    virtual void setSoleSize(const double &inLength, const double &inWidth);
    
    /**
       \brief  Get position of the ankle in the foot local coordinate frame
       
       \retval outCoordinates coordinates of the ankle joint center
    */
    virtual void getAnklePositionInLocalFrame(vector3d& outCoordinates) const;

    /**
       \brief  Set position of the ankle in the foot local coordinate frame

       \param inCoordinates coordinates of the ankle joint center
    */
    virtual void setAnklePositionInLocalFrame(const vector3d& inCoordinates);
    
    /**
       \brief Get position of the sole center in foot local frame of the foot

       \retval outCoordinates coordinates of the center C of the sole (see Figure) 
    */
    virtual void getSoleCenterInLocalFrame(vector3d& outCoordinates) const;

    /**
       \brief Set position of the sole center in foot local frame of the foot

       \param inCoordinates coordinates of the center C of the sole 
       (see Figure) 
    */
    virtual void setSoleCenterInLocalFrame(const vector3d& inCoordinates);

    /**
       \brief Get position of the projection of center of local frame in sole plane

       \retval outCoordinates coordinates of the projection H of the center of the local frame in the sole plane (see Figure) 
    */
    virtual void 
      getProjectionCenterLocalFrameInSole(vector3d& outCoordinates) const;


    /**
       \brief Set position of projection of center of local frame in sole plane

       \param inCoordinates coordinates of the projection H of the center 
       of the local frame in the sole plane (see Figure) 
    */
    virtual void 
      setProjectionCenterLocalFrameInSole(const vector3d& inCoordinates);

  private:
    /*! Store the ankle joint. */
    const CjrlJoint * m_Ankle;
    
    /*! Store sole size. */
    double m_SoleLength, m_SoleWidth;

    /*! Store ankle position in foot frame. */
    vector3d m_AnklePositionInFootFrame;

    /*! Store center position in foot frame. */
    vector3d m_CenterInFootFrame;

    /*! Store projection of center in sole frame. */
    vector3d m_ProjectionCenterInSoleFrame;
  };
   
};

#endif /* _DYN_JRL_JAPAN_FOOT_H_ */