/*
 * Copyright 2009, 2010,
 *
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
/* @doc Object used to handle specificities of humanoid robots
*/
#include "Debug.h"

#include "HumanoidSpecificities.h"
#include "fileReader.h"

using namespace dynamicsJRLJapan;

HumanoidSpecificities::HumanoidSpecificities()
{

  m_FootDepth[0] = m_FootDepth[1] = 0.14;
  m_FootHeight[0] = 0.137;
  m_FootHeight[1] = 0.137;
  m_FootWidth[0] = 0.24;
  m_FootWidth[1] = 0.24;
  m_UpperBodyJointNb = 0;
}


HumanoidSpecificities::~HumanoidSpecificities()
{
}

int HumanoidSpecificities::ReadXML(string &aFileName)
{
  FILE *fp;
  fp = fopen(aFileName.c_str(),"r");

  if (fp==0)
    {
      cerr << "Unable to read " << aFileName << endl;
      return -1;
    }

  char Side[2][80] = {"Right","Left"};
  if (look_for(fp,"Humanoid"))
    {
      ODEBUG("Found Humanoid");
      if (look_for(fp,"Feet"))
	{
	  ODEBUG("Found Feet");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{
		  ODEBUG("Found Feet Side: " << Side[i]);
		  if (look_for(fp,"SizeX"))
		    {
		      fscanfd(fp,&m_FootDepth[i]);
		      ODEBUG("Found SizeX: " << m_FootDepth[i]);
		    }

		  if (look_for(fp,"SizeY"))
		    {
		      fscanfd(fp,&m_FootWidth[i]);
		      ODEBUG("Found SizeY: " << m_FootWidth[i]);
		    }

		  if (look_for(fp,"SizeZ"))
		    {
		      fscanfd(fp,&m_FootHeight[i]);
		      ODEBUG("Found SizeZ: " << m_FootHeight[i]);
		    }

		  if (look_for(fp,"AnklePosition"))
		    {
		      fscanfd(fp, &m_AnklePosition[i][0]);
		      fscanfd(fp, &m_AnklePosition[i][1]);
		      fscanfd(fp, &m_AnklePosition[i][2]);
		      ODEBUG("AnklePos: "
			      << m_AnklePosition[i][0] << " "
			      << m_AnklePosition[i][1] << " "
			      << m_AnklePosition[i][2]);
		    }
		  if (look_for(fp,"JointNb"))
		    {
		      fscanfi(fp, &m_FeetJointNb[i]);
		      ODEBUG("JointNb: " << m_FeetJointNb[i]);
		    }
		  if (look_for(fp,"JointsID"))
		    {
		      int aJoint;
		      for(int j=0;j<m_FeetJointNb[i];j++)
			{
			  fscanfi(fp,&aJoint);
			  m_FeetJoints[i].insert( m_FeetJoints[i].end(),aJoint);
			}
		      ODEBUG("Joints :");
		      for(unsigned int j=0;j<m_FeetJoints[i].size();j++)
			{
			  ODEBUG(m_FeetJoints[i][j]);
			}

		    }

		}
	    }
	}
      else LTHROW("No feet in humanoid specificities file.");

      if (look_for(fp,"Waist"))
	{
	  ODEBUG("Waist");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{
		  ODEBUG(Side[i]);
		  if (look_for(fp,"WaistToHip"))
		    {
		      fscanfd(fp,&m_WaistToHip[i][0]);
		      fscanfd(fp,&m_WaistToHip[i][1]);
		      fscanfd(fp,&m_WaistToHip[i][2]);
		    }
		}
	    }
	  if (look_for(fp,"JointNb"))
	    {
	      fscanfi(fp, &m_WaistJointNb);
	      ODEBUG("JointNb: " << m_WaistJointNb);
	    }
	  if (look_for(fp,"JointsID"))
	    {
	      int aJoint;
	      for(int j=0;j<m_WaistJointNb;j++)
		{
		  fscanfi(fp,&aJoint);
		  m_WaistJoints.insert( m_WaistJoints.end(),aJoint);
		}
	      ODEBUG("Joints :");
	      for(unsigned int j=0;j<m_WaistJoints.size();j++)
		{
		  ODEBUG(m_WaistJoints[j]);
		}

	    }

	}
      if (look_for(fp,"Legs"))
	{
	  ODEBUG("Legs");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{
		  ODEBUG(Side[i]);

		  if (look_for(fp,"HipLength"))
		    {
		      fscanfd(fp,&m_HipLength[i][0]);
		      fscanfd(fp,&m_HipLength[i][1]);
		      fscanfd(fp,&m_HipLength[i][2]);
		      ODEBUG("Found HipLength: " << m_HipLength[i]);
		    }

		  if (look_for(fp,"FemurLength"))
		    {
		      fscanfd(fp,&m_FemurLength[i]);
		      ODEBUG("Found FemurLength: " << m_FemurLength[i]);
		    }

		  if (look_for(fp,"TibiaLength"))
		    {
		      fscanfd(fp,&m_TibiaLength[i]);
		      ODEBUG("Found TibiaLength: " << m_TibiaLength[i]);
		    }

		  if (look_for(fp,"JointNb"))
		    {
		      fscanfi(fp, &m_LegsJointNb[i]);
		      ODEBUG("JointNb: " << m_LegsJointNb[i]);
		    }
		  if (look_for(fp,"JointsID"))
		    {
		      int aJoint;
		      for(int j=0;j<m_LegsJointNb[i];j++)
			{
			  fscanfi(fp,&aJoint);
			  m_LegsJoints[i].insert( m_LegsJoints[i].end(),aJoint);
			}
		      ODEBUG("Joints :");
		      for(unsigned int j=0;j<m_LegsJoints[i].size();j++)
			{
			  ODEBUG(m_LegsJoints[i][j]);
			}

		    }
		}

	    }
	}
      if (look_for(fp,"Hands"))
	{
	  ODEBUG("Hands");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{
		  if (look_for(fp,"Center"))
		    {
		      fscanfd(fp,&m_Hands.Center[i][0]);
		      fscanfd(fp,&m_Hands.Center[i][1]);
		      fscanfd(fp,&m_Hands.Center[i][2]);
		      ODEBUG("Found Hands Center: " << m_Hands.Center[i]);
		    }

		  if (look_for(fp,"okayAxis"))
		    {
		      fscanfd(fp,&m_Hands.okayAxis[i][0]);
		      fscanfd(fp,&m_Hands.okayAxis[i][1]);
		      fscanfd(fp,&m_Hands.okayAxis[i][2]);

		      ODEBUG("Found m_Hands.okayAxis: " << m_Hands.okayAxis[i]);
		    }

		  if (look_for(fp,"showingAxis"))
		    {
		      fscanfd(fp,&m_Hands.showingAxis[i][0]);
		      fscanfd(fp,&m_Hands.showingAxis[i][1]);
		      fscanfd(fp,&m_Hands.showingAxis[i][2]);

		      ODEBUG("Found m_Hands.showingAxis: " << m_Hands.showingAxis[i]);
		    }
		  if (look_for(fp,"palmAxis"))
		    {

		      fscanfd(fp,&m_Hands.palmAxis[i][0]);
		      fscanfd(fp,&m_Hands.palmAxis[i][1]);
		      fscanfd(fp,&m_Hands.palmAxis[i][2]);

		      ODEBUG("JointNb: " << m_Hands.palmAxis[i]);
		    }
		}
	    }

	}

      if (look_for(fp,"Wrists"))
	{
	  ODEBUG("Wrists");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{

		  if (look_for(fp,"JointNb"))
		    {
		      fscanfi(fp, &m_WristsJointNb[i]);
		      ODEBUG("JointNb: " << m_WristsJointNb[i]);
		    }

		  if (look_for(fp,"JointsID"))
		    {
		      int aJoint;
		      for(int j=0;j<m_WristsJointNb[i];j++)
			{
			  fscanfi(fp,&aJoint);
			  m_WristsJoints[i].insert( m_WristsJoints[i].end(),aJoint);
			}
		      ODEBUG("Joints :");
		      for(unsigned int j=0;j<m_WristsJoints[i].size();j++)
			{
			  ODEBUG(m_WristsJoints[i][j]);
			}

		    }

		}
	    }
	}

      if (look_for(fp,"Arms"))
	{
	  ODEBUG("Arms");
	  for(int i=0;i<2;i++)
	    {
	      if (look_for(fp,Side[i]))
		{
		  if (look_for(fp,"UpperArmLength"))
		    {
		      fscanfd(fp,&m_UpperArmLength[i]);
		      ODEBUG("Found UpperArmLength: " << m_UpperArmLength[i]);
		    }

		  if (look_for(fp,"ForeArmLength"))
		    {
		      fscanfd(fp,&m_ForeArmLength[i]);
		      ODEBUG("Found ForeArmLength: " << m_ForeArmLength[i]);
		    }

		  if (look_for(fp,"JointNb"))
		    {
		      fscanfi(fp, &m_ArmsJointNb[i]);
		      ODEBUG("JointNb: " << m_ArmsJointNb[i]);
		    }
		  if (look_for(fp,"JointsID"))
		    {
		      int aJoint;
		      for(int j=0;j<m_ArmsJointNb[i];j++)
			{
			  fscanfi(fp,&aJoint);
			  m_ArmsJoints[i].insert( m_ArmsJoints[i].end(),aJoint);
			}
		      ODEBUG("Joints :");
		      for(unsigned int j=0;j<m_ArmsJoints[i].size();j++)
			{
			  ODEBUG(m_ArmsJoints[i][j]);
			}

		    }
		}
	    }
	}


      // Look for the head information.
      if (look_for(fp,"Head"))
	{

	  if (look_for(fp,"JointNb"))
	    {
	      fscanfi(fp, &m_HeadJointNb);
	      ODEBUG("JointNb: " << m_HeadJointNb);
	    }
	  if (look_for(fp,"JointsID"))
	    {
	      int aJoint;
	      for(int j=0;j<m_HeadJointNb;j++)
		{
		  fscanfi(fp,&aJoint);
		  m_HeadJoints.insert( m_HeadJoints.end(),aJoint);
		}
	      ODEBUG("Joints :");
	      for(unsigned int j=0;j<m_HeadJoints.size();j++)
		{
		  ODEBUG(m_HeadJoints[j]);
		}

	    }
	}

      // Look for the head information.
      if (look_for(fp,"Chest"))
	{

	  if (look_for(fp,"JointNb"))
	    {
	      fscanfi(fp, &m_ChestJointNb);
	      ODEBUG("JointNb: " << m_ChestJointNb);
	    }
	  if (look_for(fp,"JointsID"))
	    {
	      int aJoint;
	      for(int j=0;j<m_ChestJointNb;j++)
		{
		  fscanfi(fp,&aJoint);
		  m_ChestJoints.insert( m_ChestJoints.end(),aJoint);
		}
	      ODEBUG("Joints :");
	      for(unsigned int j=0;j<m_ChestJoints.size();j++)
		{
		  ODEBUG(m_ChestJoints[j]);
		}

	    }
	}
    }
  else LTHROW("Did not find any humanoid in humanoid specificities file.");
  fclose(fp);
  return 0;
}

void HumanoidSpecificities::Display()
{

  string Side[2] = {"Right" , "Left"};

  for(int i=0;i<2;i++)
    {
      cout << "Size of the " << Side[i] <<
	" foot:" << endl;

      cout << "Width: " <<m_FootWidth[i] << " Height: " << m_FootHeight[i] << endl;

      cout << "Ankle position of foot " << Side[i]  << endl;

      for(int j=0;j<3;j++)
	cout << m_AnklePosition[i][j] << " ";
      cout << endl;

      cout << "Number of joints :" << m_FeetJointNb[i] << endl;
      for(int j=0;j<m_FeetJointNb[i];j++)
	cout << m_FeetJoints[i][j] << " ";
      cout << endl;


      cout << "Size of the " << Side[i]
	   << " Tibia: " << m_TibiaLength[i] << endl;
      cout << "Size of the " << Side[i]
	   << " Femur: " << m_FemurLength[i] << endl;

      cout << "Number of joints :" << m_LegsJointNb[i] << endl;
      for(int j=0;j<m_LegsJointNb[i];j++)
	cout << m_LegsJoints[i][j] << " ";
      cout << endl;

      cout << "Size of the " << Side[i]
	   << " UpperArm: " << m_UpperArmLength[i] << endl;
      cout << "Size of the " << Side[i]
	   << " ForeArm: " << m_ForeArmLength[i] << endl;

      cout << "Number of joints :" << m_ArmsJointNb[i] << endl;
      for(int j=0;j<m_ArmsJointNb[i];j++)
	cout << m_ArmsJoints[i][j] << " ";
      cout << endl;

      cout << "Size of the foot " << Side[i]
	   << " Width: " << m_FootWidth[i]
	   << " Height: " << m_FootHeight[i] << endl;

      cout << "Position of the hip according to the waist  " << Side[i] << " " ;
      for(int j=0;j<3;j++)
	cout << m_WaistToHip[i][j] << " ";
      cout << endl;

    }


}

int HumanoidSpecificities::GetFootSize(int WhichFoot, double &Depth,
				       double & Width, double &Height)
{
  Width = 0.0;
  Height = 0.0;
  Depth = 0.0;

  if (WhichFoot==-1)
    {
      Depth = m_FootDepth[0];
      Width = m_FootWidth[0];
      Height = m_FootHeight[0];
    }
  else if (WhichFoot==1)
    {
      Depth = m_FootDepth[1];
      Width = m_FootWidth[1];
      Height = m_FootHeight[1];
    }
  else
    return -1;

  return 0;
}


double HumanoidSpecificities::GetTibiaLength(int WhichSide)
{
  if (WhichSide==-1)
    return m_TibiaLength[0];
  return m_TibiaLength[1];
}

double HumanoidSpecificities::GetFemurLength(int WhichSide)
{
  if (WhichSide==-1)
    return m_FemurLength[0];
  return m_FemurLength[1];
}

double HumanoidSpecificities::GetUpperArmLength(int WhichSide)
{
  if (WhichSide==-1)
    return m_UpperArmLength[0];
  return m_UpperArmLength[1];
}

double HumanoidSpecificities::GetForeArmLength(int WhichSide)
{
  if (WhichSide==-1)
    return m_ForeArmLength[0];
  return m_ForeArmLength[1];
}

void HumanoidSpecificities::GetAnklePosition(int WhichSide,double AnklePosition[3])
{
  int r=1;
  if (WhichSide==-1)
    r=0;
  for (int i=0;i<3;i++)
    AnklePosition[i] = m_AnklePosition[r][i];
}

void HumanoidSpecificities::GetWaistToHip(int WhichSide,double WaistToHip[3])
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  for(int i=0;i<3;i++)
    WaistToHip[i] = m_WaistToHip[r][i];
}

void HumanoidSpecificities::GetHipLength(int WhichSide,double HipLength[3] )
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  for(int i=0;i<3;i++)
    HipLength[i] = m_HipLength[r][i];

}

int HumanoidSpecificities::GetArmJointNb(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_ArmsJointNb[r];
}

const std::vector<int> & HumanoidSpecificities::GetArmJoints(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_ArmsJoints[r];
}

const std::vector<int> &HumanoidSpecificities::GetWrists(int WhichSide)
{
  int r=1;

  if (WhichSide==-1)
    r=0;

  return m_WristsJoints[r];
}

int HumanoidSpecificities::GetLegJointNb(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_LegsJointNb[r];
}

const std::vector<int> & HumanoidSpecificities::GetLegJoints(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_LegsJoints[r];
}

int HumanoidSpecificities::GetFootJointNb(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_FeetJointNb[r];
}

const std::vector<int> & HumanoidSpecificities::GetFootJoints(int WhichSide)
{
  int r=1;
  if (WhichSide==-1)
    r=0;

  return m_FeetJoints[r];
}

int HumanoidSpecificities::GetHeadJointNb()
{
  return m_HeadJointNb;
}

const std::vector<int> & HumanoidSpecificities::GetHeadJoints()
{
  return m_HeadJoints;
}

int HumanoidSpecificities::GetChestJointNb()
{
  return m_ChestJointNb;
}

const std::vector<int> & HumanoidSpecificities::GetChestJoints()
{
  return m_ChestJoints;
}

int HumanoidSpecificities::InitUpperBodyJoints()
{
  m_UpperBodyJointNb = m_HeadJointNb + m_ChestJointNb +
    m_ArmsJointNb[0] + m_ArmsJointNb[1];
  int lindex =0;

  m_UpperBodyJoints.resize(m_UpperBodyJointNb);

  for(int i=0;i<m_HeadJointNb;i++)
    m_UpperBodyJoints[lindex++] = m_HeadJoints[i];

  for(int i=0;i<m_ChestJointNb;i++)
    m_UpperBodyJoints[lindex++] = m_ChestJoints[i];

  for(int i=0;i<m_ArmsJointNb[0];i++)
    m_UpperBodyJoints[lindex++] = m_ArmsJoints[0][i];

  for(int i=0;i<m_ArmsJointNb[1];i++)
    m_UpperBodyJoints[lindex++] = m_ArmsJoints[1][i];

  return 0;
}

int HumanoidSpecificities::GetUpperBodyJointNb()
{
  if (m_UpperBodyJointNb==0)
    {
      InitUpperBodyJoints();
    }
  return m_UpperBodyJointNb;
}


const std::vector<int> & HumanoidSpecificities::GetUpperBodyJoints()
{
  if (m_UpperBodyJointNb==0)
    {
      InitUpperBodyJoints();
    }
  return m_UpperBodyJoints;
}

int HumanoidSpecificities::GetWaistJointNb()
{
  return m_WaistJointNb;
}

const std::vector<int> & HumanoidSpecificities::GetWaistJoints()
{
  return m_WaistJoints;
}

const HandsData & HumanoidSpecificities::GetHandsData()
{
  return m_Hands;
}
