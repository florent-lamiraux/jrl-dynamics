/* @doc \file Object to generate a file following VRML 1.0 format.

   Copyright (c) 2010
   @author Olivier Stasse
   JRL-Japan, CNRS/AIST
 
   All rights reserved.

   Please refers to file License.txt for details on the license.

*/
/* System includes */
#include <iostream>
#include <fstream>

/*! Includes of the framework */
#include "GenerateRobotForVRML2.h"

using namespace std;

namespace dynamicsJRLJapan {

  Tools::GenerateRobotForVRML2::GenerateRobotForVRML2()
  {
  }


  Tools::GenerateRobotForVRML2::~GenerateRobotForVRML2()
  {
  }


  void Tools::GenerateRobotForVRML2::AxisAngle(matrix4d &data,
					       vector3d &axis,
					       double &angle) const
  {
    MAL_S3_VECTOR_FILL(axis,0.0);
    angle = 0;
    double x,y,z;
    //    cout << "data:" << data << endl;

    double  d0 = MAL_S4x4_MATRIX_ACCESS_I_J(data,0,0);
    double d01 = MAL_S4x4_MATRIX_ACCESS_I_J(data,0,1);
    double d02 = MAL_S4x4_MATRIX_ACCESS_I_J(data,0,2);
    double d10 = MAL_S4x4_MATRIX_ACCESS_I_J(data,1,0);
    double  d1 = MAL_S4x4_MATRIX_ACCESS_I_J(data,1,1);
    double d12 = MAL_S4x4_MATRIX_ACCESS_I_J(data,1,2);
    double d20 = MAL_S4x4_MATRIX_ACCESS_I_J(data,2,0);
    double d21 = MAL_S4x4_MATRIX_ACCESS_I_J(data,2,1);
    double  d2 = MAL_S4x4_MATRIX_ACCESS_I_J(data,2,2);

    double epsilon = 0.01; // margin to allow for rounding errors
    double epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees
    
    if ((fabs(d01 - d10)<epsilon) &&
	(fabs(d02 - d20)<epsilon) &&
	(fabs(d12 - d21)<epsilon))
      {
	if ((fabs(d01 - d10)<epsilon2) &&
	    (fabs(d02 - d20)<epsilon2) &&
	    (fabs(d12 - d21)<epsilon2))
	  {
	    MAL_S3_VECTOR_ACCESS(axis,0)= 0.0;MAL_S3_VECTOR_ACCESS(axis,1)= 1.0;MAL_S3_VECTOR_ACCESS(axis,2)= 0.0;
	    angle=0.0;
	    return;
	  }
	angle = M_PI;
	double xx = (d0+1)/2;
	double yy = (d1+1)/2;
	double zz = (d2+1)/2;
	double xy = (d01+d10)/4;
	double xz = (d02+d20)/4;
	double yz = (d12+d21)/4;
	if ((xx > yy) && (xx > zz)) { // m[0][0] is the largest diagonal term
	  if (xx< epsilon) {
	    x = 0;
	    y = 0.7071;
	    z = 0.7071;
	  } else {
	    x = sqrt(xx);
	    y = xy/x;
	    z = xz/x;
	  }
	} else if (yy > zz) { // m[1][1] is the largest diagonal term
	  if (yy< epsilon) {
	    x = 0.7071;
	    y = 0;
	    z = 0.7071;
	  } else {
	    y = sqrt(yy);
	    x = xy/y;
	    z = yz/y;
	  }	
	} else { // m[2][2] is the largest diagonal term so base result on this
	  if (zz< epsilon) {
	    x = 0.7071;
	    y = 0.7071;
	    z = 0;
	  } else {
	    z = sqrt(zz);
	    x = xz/z;
	    y = yz/z;
	  }
	}
	MAL_S3_VECTOR_ACCESS(axis,0) = x;
	MAL_S3_VECTOR_ACCESS(axis,1) = y;
	MAL_S3_VECTOR_ACCESS(axis,2) = z;
	return;
      }
    double s = sqrt( (d21 - d12)*(d21 -d12 ) + 
		     (d02 - d20)*(d02 - d20) +
		     (d10 - d01)*(d10 - d01));
    if (fabs(s)<0.001) s =1.0;
    double r = (d0+d1+d2-1.0)*0.5 ;
    //    cout << "r:" << r << endl;
    angle = acos(r);
    //angle = 0;
    //    cout << "angle:" << angle<<endl;
    MAL_S3_VECTOR_ACCESS(axis,0) = (d21 - d12)/s;
    MAL_S3_VECTOR_ACCESS(axis,1) = (d02 - d20)/s;
    MAL_S3_VECTOR_ACCESS(axis,2) = (d10 - d01)/s;

    return;
  }

  /*! \brief Generate hard coded robot */
  void Tools::GenerateRobotForVRML2::GenerateRobot(std::string &RobotName,
			    CjrlHumanoidDynamicRobot *aHDR)
  {
    m_HDR = aHDR;
    ofstream aof;
    string FileName = RobotName + "_v1.wrl";
    aof.open((char *)FileName.c_str(),ofstream::out);

    if(!aof.is_open())
      return;

    aof << "#VRML V2.0 utf8" << endl;
    aof << endl;
    aof << "Separator {" << endl;
    string shifttab="  ";
    GenerateBodies(aof,shifttab);
    aof << "}" << endl;
    aof.close();
  }


  void Tools::GenerateRobotForVRML2::GenerateBodies(ostream &os,
						    string shifttab)
  {
    CjrlJoint *RootJoint = m_HDR->rootJoint();
    string lshifttab = shifttab+"  ";
    unsigned int gindex=0;

    GenerateJoint(RootJoint,os,shifttab,gindex);
    os << shifttab << "children [" << endl;
    GenerateBody(RootJoint,os, lshifttab,gindex);
    os << shifttab << "]" << endl;

  }
  
  /*! \brief Take the links towards the geometrical information */
  void Tools::GenerateRobotForVRML2::SetAccessToData(std::vector<std::string> &AccessToData)
  {
    m_AccessToData = AccessToData;
  }

  void Tools::GenerateRobotForVRML2::GenerateJoint(CjrlJoint *aJoint, 
					     ostream &os,
					     string shifttab, unsigned int &gindex)
  {

    os << shifttab << "Transform { "<< endl;
    
    MAL_S4x4_MATRIX(aTransformation,double);
    aTransformation = aJoint->currentTransformation();
    os << shifttab << "  translation ";
    for(unsigned int i=0;i<3;i++)
      os << MAL_S4x4_MATRIX_ACCESS_I_J(aTransformation,i,3) << " ";
    os << endl;
    os << shifttab << "  rotation " ;
    vector3d laxis;
    double angle;
    AxisAngle(aTransformation,laxis,angle);
    for(unsigned int i=0;i<3;i++)
      os << MAL_S3_VECTOR_ACCESS(laxis,i) << " ";
    
    os << angle << endl;

  }  


  void Tools::GenerateRobotForVRML2::SetPathToModelFiles(std::string &Path)
  {
    m_PathToModelFiles = Path;
  }
  void Tools::GenerateRobotForVRML2::CopyGeometricInformation(ostream &os,
							      string FileName)
  {
    string FinalFileName = m_PathToModelFiles+ FileName;
    std::ifstream geometricFile((char *)FinalFileName.c_str(),ifstream::in);
    if (!geometricFile.is_open())
      {
	cerr << "Unable to open " << FinalFileName << endl;
	return;
      }
    // Go through the geometric file
    while(!geometricFile.eof())
      {
	char geometricline[25536];
	geometricFile.getline(geometricline,25536);
	string aline(geometricline);
	os.write(geometricline,aline.length());
      }
    geometricFile.close();

  }

  void Tools::GenerateRobotForVRML2::GenerateBody(CjrlJoint *aJoint, 
						  ostream &os,
						  string shifttab,
						  unsigned int &gindex)
  {

    CjrlBody *aBody= aJoint->linkedBody();
    if (aBody==0)
      return;

    GenerateJoint(aJoint,os,shifttab,gindex);

        // Geometric file.
    os << shifttab << "  children [" << endl
       << shifttab << "    Inline {"<< endl
       << shifttab << "      url \"" 
       << m_AccessToData[gindex] << "\"" << endl
       << shifttab << "    }" << endl;
    
      //    CopyGeometricInformation(os,m_AccessToData[gindex]);
    os << shifttab << "  ]"<< endl;;
    os << shifttab << "}" << endl;  

    gindex++;
    // Call the sons.
    for(unsigned int i=0;i<aJoint->countChildJoints();i++)
      {
	GenerateBody(aJoint->childJoint(i),os,shifttab,gindex);
      }

  }


}
