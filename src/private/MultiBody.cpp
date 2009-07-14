/* Multibody object used to compute :
   - Center Of Mass,
   - Zero Momentum Point.

   OS: Updated the names of the contributors, the documentation
   and added a sample file for WalkPlugin
   OS (21/12/2006): removed any reference to non-homogeneous matrix
   library.
   OS (10/01/2007): Put the abstract layer for small matrix library.
   
   Copyright (c) 2005-2006, 
   Adrien Escande,
   Abderrahmane Kheddar,  
   Olivier Stasse,
   Ramzi Sellouati
   
   JRL-Japan, CNRS/AIST

   Please refers to file License.txt for details on the license.

*/
#include "MultiBody.h"
#include "SpiritVRMLReader.h"

using namespace dynamicsJRLJapan;

// surcharge operateur
bool dynamicsJRLJapan::operator==(const dynamicsJRLJapan::appariement a1, 
				        const dynamicsJRLJapan::appariement a2) 
{
  return (a1.corps==a2.corps);
}



// Caution: This operator is specific to OpenGL matrices: it transposes the
// matrix before multiplication.
MAL_S3_VECTOR(,double) operator * (double* m, MAL_S3_VECTOR(,double) v) 
{
  MAL_S3_VECTOR(,double) result;

  result[0] = m[0]*v[0] + m[4]*v[1] +  m[8]*v[2]+ m[12];
  result[1] = m[1]*v[0] + m[5]*v[1] +  m[9]*v[2]+ m[13];
  result[2] = m[2]*v[0] + m[6]*v[1] + m[10]*v[2]+ m[14];
  return result;
}


//Pour matrice OpenGL
void Matrix2AxeAngle(float R[16],double Axis[3], double & Angle)
{
  double q[4];
  double sum_x, sum_y, sum_z, sum_w, sum_max, S;
  
  sum_w = 1 + R[0] + R[5] + R[10];
  sum_x = 1 + R[0] - R[5] - R[10];
  sum_y = 1 + R[5] - R[0] - R[10];
  sum_z = 1 + R[10] - R[0] - R[5];
  sum_max = max(max(sum_w,sum_x), max(sum_y, sum_z));
  
  if (sum_max == sum_w)
    {
      S = sqrt(sum_w)*2;
      q[0] = (R[9] - R[6])/S;
      q[1] = (R[2] - R[8])/S;
      q[2] = (R[4] - R[1])/S;
      q[3] = 0.25*S;
    }
  else if (sum_max == sum_x)
    {
      S  = sqrt(sum_x) * 2;
      q[0] = 0.25 * S;
      q[1] = (R[4] + R[1] ) / S;
      q[2] = (R[2] + R[8] ) / S;
      q[3] = (R[9] - R[6] ) / S;
    } 
  else if (sum_max == sum_y) 
    { 
      S  = sqrt(sum_y) * 2;
      q[0] = (R[4] + R[1] ) / S;
      q[1] = 0.25 * S;
      q[2] = (R[9] + R[6] ) / S;
      q[3] = (R[2] - R[8] ) / S;
    } 
  else
    {
      S  = sqrt(sum_z) * 2;
      q[0] = (R[2] + R[8] ) / S;
      q[1] = (R[9] + R[6] ) / S;
      q[2] = 0.25 * S;
      q[3] = (R[4] - R[1] ) / S;
    } 

  //passage des quaternions aux axes et angles
  double sum;
  double cos_a,  sin_a;

  // Normalize
  sum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
  sum = sqrt(sum);

  if (sum==0)
    {
      Angle = Axis[0] = Axis[1] = Axis[2] = 0.0;
      return;
    }

  for(int i=0;i<4;i++)
    q[i] /= sum;

  cos_a = q[3];
  Angle = acos( cos_a ) * 2;
  sin_a = sqrt( 1.0 - cos_a * cos_a );
  if ( fabs( sin_a ) < 0.0005 ) 
    sin_a = 1;

  Axis[0] = -q[0] / sin_a;
  Axis[1] = -q[1] / sin_a;
  Axis[2] = -q[2] / sin_a;
}



MultiBody::MultiBody(void)
{
    masse = 0.0;
}


MultiBody::~MultiBody(void)
{
  for(unsigned int li=0;li<listeLiaisons.size();li++)
    delete listeLiaisons[li].aJoint;

  for(unsigned int li=0;li<listeCorps.size();li++)
    delete listeCorps[li];
}

void MultiBody::ajouterCorps(Body &b)
{
  masse += b.getMasse();
  listeCorps.push_back(&b);
  liaisons.push_back(vector<appariement>());
}

Body * MultiBody::dernierCorps()
{
  return listeCorps[listeCorps.size()-1];
}

void MultiBody::ajouterLiaison(Body &corps1, Body &corps2, internalLink & l)
{
  //recherche de l'index du premier corps dans listeCorps
  unsigned int index1;
  for (index1 = 0; index1<listeCorps.size(); index1++) {
    if (corps1.getLabel() == listeCorps[index1]->getLabel()) {
      break;
    }
  }

  //recherche de l'index du deuxieme corps dans listeCorps
  unsigned int index2;
  for (index2 = 0; index2<listeCorps.size(); index2++) {
    if (corps2.getLabel() == listeCorps[index2]->getLabel()) {
      break;
    }
  }

  l.indexCorps1 = index1;
  l.indexCorps2 = index2;
  listeLiaisons.push_back(l);		//ajout de la liaison a la liste

  //creation de l'appariement corps2, liaison
  appariement a2 = {index2, listeLiaisons.size()-1};
  liaisons[index1].push_back(a2);	//et ajout
  //creation de l'appariement corps1, liaison
  appariement a1 = {index1, listeLiaisons.size()-1};
  liaisons[index2].push_back(a1);	//et ajout
}

void MultiBody::ajouterLiaisonFixe(Body &corps1, Body &corps2, 
				   MAL_S3_VECTOR(,double) translationStat, 
				   MAL_S3_VECTOR(,double) axeRotationStat, 
				   double angleRotationStat)
{
  //recherche de l'index du premier corps dans listeCorps
  unsigned int index1;
  for (index1 = 0; index1<listeCorps.size(); index1++) {
    if (corps1.getLabel() == listeCorps[index1]->getLabel()) {
      break;
    }
  }

  //recherche de l'index du deuxieme corps dans listeCorps
  unsigned int index2;
  for (index2 = 0; index2<listeCorps.size(); index2++) {
    if (corps2.getLabel() == listeCorps[index2]->getLabel()) {
      break;
    }
  }

  //creation de la liaison
  JointPrivate * aJoint= new JointPrivate(JointPrivate::FIX_JOINT,axeRotationStat,angleRotationStat,translationStat);

  internalLink l = {cptLiaison++,aJoint,index1, index2};
  listeLiaisons.push_back(l);		//ajout de la liaison a la liste
  //creation de l'appariement corps2, liaison
  appariement a2 = {index2, listeLiaisons.size()-1};
  liaisons[index1].push_back(a2);	//et ajout
  //creation de l'appariement corps1, liaison
  appariement a1 = {index1, listeLiaisons.size()-1};
  liaisons[index2].push_back(a1);	//et ajout
}

void MultiBody::supprimerLiaisonEntre(Body &corps1, Body &corps2)
{
  //recherche de l'index du premier corps dans listeCorps
  unsigned int c1;
  for (c1 = 0; c1<listeCorps.size(); c1++) {
    if (corps1.getLabel() == listeCorps[c1]->getLabel()) {
      break;
    }
  }

  //recherche de l'index du deuxieme corps dans listeCorps
  unsigned int c2;
  for (c2 = 0; c2<listeCorps.size(); c2++) {
    if (corps2.getLabel() == listeCorps[c2]->getLabel()) {
      break;
    }
  }

  appariement a1 = {c2, 0};
  appariement a2 = {c1, 0};
  //recherche de l'iterateur de la liaison parmis les liaisons de corps1
  vector<appariement>::iterator it = find(liaisons[c1].begin(), liaisons[c1].end(), a1);
  int index = it->liaison;
  listeLiaisons[index].indexCorps1 = -1;
  listeLiaisons[index].indexCorps2 = -1;
  listeLiaisons[index].label = -1;
  liaisons[c1].erase(it);		//et on l'enleve
  //recherche de l'iterateur de la liaison parmis les liaisons de corps2
  it = find(liaisons[c2].begin(), liaisons[c2].end(), a2);
  liaisons[c2].erase(it);		//et on l'enleve
}

void MultiBody::supprimerLiaison(int index)
{
  int c1 = listeLiaisons[index].indexCorps1;
  int c2 = listeLiaisons[index].indexCorps2;
  appariement a1 = {c2, 0};
  appariement a2 = {c1, 0};
  //recherche de l'iterateur de la liaison parmis les liaisons de corps1
  vector<appariement>::iterator it = find(liaisons[c1].begin(), liaisons[c1].end(), a1);
  listeLiaisons[index].indexCorps1 = -1;
  listeLiaisons[index].indexCorps2 = -1;
  listeLiaisons[index].label = -1;
  liaisons[c1].erase(it);		//et on l'enleve
  //recherche de l'iterateur de la liaison parmis les liaisons de corps2
  it = find(liaisons[c2].begin(), liaisons[c2].end(), a2);
  liaisons[c2].erase(it);		//et on l'enleve
}

void MultiBody::supprimerLiaisonLabel(int label)
{
  //recherche de l'index de la liaison
  unsigned int i;
  for (i=0; i<listeLiaisons.size(); i++) {
    if (listeLiaisons[i].label == label)
      break;
  }
  if (i==listeLiaisons.size())
    return;

  supprimerLiaison(i);
}

void MultiBody::supprimerCorps(Body &b)
{
  unsigned int i;
  for (i=0; i<listeCorps.size(); i++) {
    if (listeCorps[i]->getLabel() == b.getLabel())
      break;
  }
  if (i==listeCorps.size())
    return;

  supprimerCorps(i);
}
void MultiBody::supprimerCorps(int index)
{
  int j = 0;
  for (unsigned int i=0; i<liaisons[index].size(); j++) {
    supprimerLiaison(liaisons[index][i].liaison);
    cout << "suppression liaison\n";
  }
  listeCorps[index]->setLabel(-1);

}
void MultiBody::supprimerCorpsLabel(int label)
{
  unsigned int i;
  for (i=0; i<listeCorps.size(); i++) {
    if (listeCorps[i]->getLabel() == label)
      break;
  }
  if (i==listeCorps.size())
    return;

  supprimerCorps(i);
}

void MultiBody::inverserLiaison(int i)
{
  cout << "inverserLiaison n'est pas correctement implemente" << endl;
  /*
    if ((unsigned int)i >= listeLiaisons.size()) {
    cout << i << " : indice hors limite du vecteur de liaison" << endl;
    return;
    }
    for (unsigned int j=0; j<listeLiaisons[i].listeJoints.size(); j++) {
    listeLiaisons[i].listeJoints[j].axe = -listeLiaisons[i].listeJoints[j].axe;
    }
  */
}



MAL_S3_VECTOR(,double) MultiBody::getPositionCoM(void)
{
  return (positionCoMPondere/masse);
}


void MultiBody::afficher()
{
  afficherCorps();
  afficherLiaisons();
}

void MultiBody::afficherCorps()
{
  for (unsigned int i=0; i<listeCorps.size(); i++) {
    cout << "corps "<< i << " : \n";
    listeCorps[i]->Display();
    for (unsigned int j=0; j<liaisons[i].size(); j++) {
      cout << "    lie a corps " << liaisons[i][j].corps << " par liaison " 
	   << liaisons[i][j].liaison << " (label " << listeLiaisons[liaisons[i][j].liaison].label <<")\n";
    }
    cout << "\n";
  }
  cout << "\n";
}

void MultiBody::afficherLiaisons(void) {
  for (unsigned int i=0; i<listeLiaisons.size(); i++) {
    cout << "Name: "<< listeLiaisons[i].aJoint->getName()
	 << " JointID in VRML " 
	 << listeLiaisons[i].aJoint->getIDinActuated() << " " ;
    cout << "liaison de type " << listeLiaisons[i].aJoint->type()
	 << "  label "<< listeLiaisons[i].label 
	 << "  liant le corps " 
	 << listeLiaisons[i].indexCorps1 
	 << " au corps " << listeLiaisons[i].indexCorps2 << "\n";
    cout << "translationStatique : " << endl;
    MAL_S3_VECTOR(,double) aStaticTranslation;
    listeLiaisons[i].aJoint->getStaticTranslation(aStaticTranslation);
    cout << aStaticTranslation << endl;
    if (listeLiaisons[i].aJoint->type() > 0) {
      cout << "    axe : " << endl;
      cout << listeLiaisons[i].aJoint->axe();
      cout << endl;
    }
  }
  cout << "\n";
}

void dynamicsJRLJapan::AxeAngle2Matrix(const vector3d &AnAxis, double aQuantity, matrix3d &R)
{
    const double c = cos(aQuantity);
    const double s = sin(aQuantity);
    const double v = 1.0-c;
    const double xv  = AnAxis[0]*AnAxis[0]*v;
    const double yv  = AnAxis[1]*AnAxis[1]*v;
    const double zv  = AnAxis[2]*AnAxis[2]*v;
    const double xyv = AnAxis[0]*AnAxis[1]*v;
    const double yzv = AnAxis[1]*AnAxis[2]*v;
    const double zxv = AnAxis[2]*AnAxis[0]*v;
    const double xs  = AnAxis[0]*s;
    const double ys  = AnAxis[1]*s;
    const double zs  = AnAxis[2]*s;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,0,0) = xv+c;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,0,1) = xyv - zs;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,0,2) = zxv + ys;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,1,0) = xyv + zs;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,1,1) = yv + c;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,1,2) = yzv - xs;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,2,0) = zxv - ys;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,2,1) = yzv + xs;
    MAL_S3x3_MATRIX_ACCESS_I_J(R,2,2) = zv + c;
}

void MultiBody::parserVRML(string path, string nom, const char* option)
{
  string nomWRML = path;
  nomWRML += nom;
  dynamicsJRLJapan::VRMLReader::ParseVRMLFile(this,nomWRML);
}

void MultiBody::parserVRML(string nomWRML, const char* option)
{
  dynamicsJRLJapan::VRMLReader::ParseVRMLFile(this,nomWRML);
}



int MultiBody::NbOfLinks() const
{
  return listeLiaisons.size();
}

int MultiBody::NbOfJoints() const
{
  return listeLiaisons.size();
}

double MultiBody::getMasse()
{
  return masse;
}