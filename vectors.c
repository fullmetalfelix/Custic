#include "defos.h"


//gets the norm of a vector v using the overlap matrix as metric!
//(the size has to be pSize)
double Norm(double *v)
{
	int i,j;
	double tmp=0.00f;
	double result=0.00f;
	
	/*printf("vector is\n");
	for(i=0;i<pSize;i++)
		printf("%f,",v[i]);
	printf("\n");*/
	
	for(j=0;j<pSize;j++)
	{
		tmp=0.00f;
		for(i=0;i<pSize;i++)
		{
			tmp+=v[j]*SMatrix[SMap[j+i*pSize]]*v[i];
		}
		result+=tmp;
		//printf("ASD! %f\n",result);
	}
	//printf("NORM!! %f\n",sqrt(result));
	//printf("&LOL! %d\n",&result);	
	//result=2.0f;
	result=sqrt(result);
	return result;
	
}

//normalize a vector of size pSize
int Normalize(double *vector)
{
	int i;
	double nrm=(Norm(vector)); //get the norm
	
	for(i=0;i<pSize;i++)
		vector[i]/=nrm;
	return 0;
	
}


//get the norm of a single vector in a matrix
double NormOffset(double *Matrix, int row)
{
  int i,j;
  double tmp=0.00f;
  double result=0.00f;
  
  for(j=0;j<pSize;j++)
    {
      tmp=0.00f;
      for(i=0;i<pSize;i++)
	{
	  tmp+=Matrix[row*pSize + j]*SMatrix[SMap[j+i*pSize]]*Matrix[row*pSize + i];
	}
      result+=tmp;
      //printf("ASD! %f\n",result);
    }
  
  result=sqrt(result);
  return result;
  
}

//normalize a vector in a matrix
int NormalizeOffset(double *Matrix, int row)
{
	int i;
	double nrm=(NormOffset(Matrix,row));
	
	//printf("Norm of %i = %f\n",row,nrm);
	
	for(i=0;i<pSize;i++)
		Matrix[row*pSize + i]/=nrm;
	return 0;
	
}

//feed in the index of the orbital and this will give a pointer to the atom
//that holds that orbital and the shell index of that orbital!
Atom *GetAtom(int OrbitalIndex,int *shellIdx)
{
	int i,j,k;
	Atom* atm;
	
	k=0;
	for(i=0;i<SysInfo.AtomsInSystem;i++)			//loop on all atoms 
	{
		for(j=0;j<AtomSystem[i].Type->nShells;j++)
		{
			if(k==OrbitalIndex)
			{
				*shellIdx = j;
				//printf("the atom is in the memory %d, with orbital %d\n",&AtomSystem[i],*shellIdx);
				return &AtomSystem[i];
			}
			k++;
		}
		
	}
	return 0;
	
}

//input the atom index and the orbital index(refered to the atom)
//and this gives the global index for the orbital
int GetOrbital(int AtomIndex, int OrbIndex)
{
  int i,j,k;

  k=0;
  for(i=0;i<=AtomIndex;i++)
    {
      for(j=0;j<AtomSystem[i].Type->nShells;j++)
	{
	  k++;	  
	  if( (j == OrbIndex) && (i == AtomIndex) )
	    return k-1;
	}
    }

  return 0;

}


//this computes the SQUARE of the distance between atom1 and atom2
double GetAtomDistance(Atom *atom1, Atom *atom2)
{
	double result=0.0f;
	
	result+=pow(((*atom1).Position.x-(*atom2).Position.x),2.0);
	result+=pow(((*atom1).Position.y-(*atom2).Position.y),2.0);
	result+=pow(((*atom1).Position.z-(*atom2).Position.z),2.0);
	//result=sqrt(result);
	
	if(result<RkcTolerance)
		result=0.0;
	
	return result;
	
}
//same as previous but starts from vectors instead of atoms!
double GetVectorDistance(Vector3 r1, Vector3 r2)
{
	double result=0.0;
	
	result+=pow(r1.x-r2.x,2.0);
	result+=pow(r1.y-r2.y,2.0);
	result+=pow(r1.z-r2.z,2.0);
	if(result<RkcTolerance)
		return 0.0;
	//result=sqrt(result);
	return result;
	
}

//returns the vector r1-r2
Vector3 GetVectorMinus(Vector3 r1, Vector3 r2)
{
	Vector3 r;
	
	r.x = r1.x-r2.x;
	r.y = r1.y-r2.y;
	r.z = r1.z-r2.z;
	//printf("vector minus 1: %f %f %f\n",r1.x,r1.y,r1.z);
	//printf("vector minus 2: %f %f %f\n",r2.x,r2.y,r2.z);
	//printf("vector minus 3: %f %f %f\n",r.x,r.y,r.z);
	return r;
}

//returns the vector r1-r2 where r1(r2) is the position of atom1 (atom2)
Vector3 GetVectorMinusA(Atom *atom1, Atom *atom2)
{
	Vector3 r1 = (*atom1).Position;
	Vector3 r2 = (*atom2).Position;
	Vector3 r = GetVectorMinus(r1, r2);
	//printf("vector minusA: %f %f %f\n",r.x,r.y,r.z);
	return r;
}

//returns the c-th component of r
double VectorComponent(Vector3 r, int c)
{
	
	//printf("vector comp: %f %f %f\n",r.x,r.y,r.z);
	
	if(c==0)
		return r.x;
	if(c==1)
		return r.y;
	if(c==2)
		return r.z;
	
	printf("ERROR!!!");
	return 0.0f;
	
}

//computes Rkc^2
double GetRkc(Atom *atom1, Atom *atom2, Atom *atomc, double a, double b)
{
	double result=0.0f;
	Vector3 Rk;
	
	Rk.x=((*atom1).Position.x*a+(*atom2).Position.x*b)/(a+b);
	Rk.y=((*atom1).Position.y*a+(*atom2).Position.y*b)/(a+b);
	Rk.z=((*atom1).Position.z*a+(*atom2).Position.z*b)/(a+b);
	
	result= pow(((*atomc).Position.x - Rk.x),2.0f);
	result+=pow(((*atomc).Position.y - Rk.y),2.0f);
	result+=pow(((*atomc).Position.z - Rk.z),2.0f);
	
	if(result<RkcTolerance)
		result=0.0f;
	
	return result;
}


//compute the square norm of vector P12 = Sum[ (ax1 + bx2)/(a+b) , xyz]
double GetR12Norm2(Atom *atom1, Atom *atom2, double a, double b)
{
	double result=0.0f;
	Vector3 Rk;
	
	Rk.x=((*atom1).Position.x*a+(*atom2).Position.x*b)/(a+b);
	Rk.y=((*atom1).Position.y*a+(*atom2).Position.y*b)/(a+b);
	Rk.z=((*atom1).Position.z*a+(*atom2).Position.z*b)/(a+b);
	
	result = Rk.x*Rk.x + Rk.y*Rk.y + Rk.z*Rk.z;
	if(result<RkcTolerance)
		result=0.0f;
	return result;
	
}

//returns a vector with the components of P12
Vector3 GetR12Vect(Atom *atom1, Atom *atom2, double a, double b)
{
	//double result=0.0f;
	Vector3 Rk;
	
	Rk.x=((*atom1).Position.x*a+(*atom2).Position.x*b)/(a+b);
	Rk.y=((*atom1).Position.y*a+(*atom2).Position.y*b)/(a+b);
	Rk.z=((*atom1).Position.z*a+(*atom2).Position.z*b)/(a+b);
	
	return Rk;
	
}



void DisplaceAtom(int lambda, int xyz, double dx)
{

  if(xyz == 0)
    AtomSystem[lambda].Position.x += dx;
  if(xyz == 1)
    AtomSystem[lambda].Position.y += dx;
  if(xyz == 2)
    AtomSystem[lambda].Position.z += dx;
}



//copy the values from matIn to matOut
int MatrixCopy(double *matIn, double *matOut, int SizeX, int SizeY)
{
	int i;
	for(i=0;i<SizeY*SizeX;i++)
	{
		matOut[i]=matIn[i];
	}
	return 0;
	
}
//same function but with integer vectors
int MatrixCopyInt(int *matIn, int *matOut, int SizeX, int SizeY)
{
	int i;
	for(i=0;i<SizeY*SizeX;i++)
	{
		matOut[i]=matIn[i];
	}
	return 0;
	
}



