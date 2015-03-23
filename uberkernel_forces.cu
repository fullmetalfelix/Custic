// -*- c++ -*-


//*** ***** FORCES EVALUATION KERNELS ***** ***//
//***********************************************

__device__ void CoulombForce(float r[], float dist, float Cc, float Force[])
{

 
  float qq = dist * dist * dist;
  qq = Cc *14.3996368420f / qq;

  #pragma unroll
  for(int c=0;c<3;c++)
    {
      //r[c] *= dist;
      Force[c] -= r[c]*qq;
    } 
  
}

__device__ void PairForce(float r[3], float dist, int i, int j, float Force[3])
{

  float4 p =  ShellParameters_d[ max(i,j)*(max(i,j)+1)/2 + min(i,j) ];  //!! CONST MEMORY READ

  //this is -grad(E)
  float qq = (-p.x*expf(-dist/p.y)/p.y + p.z*p.w*powf(dist,p.w-1)) / dist;
 
  //if(threadIdx.x == 0)
  //printf("pair interacton %f\n",qq*r[0]);
 
  #pragma unroll
  for(int c=0;c<3;c++)
    Force[c] += r[c]*qq ;
 
}



__global__ void Calc_Forces( float4 *pos, short *AtomType_d,  float3 *forces )
{
  
  __shared__ float4 spos[BLOCKSIZE];           //shared positions
  __shared__ short  types[BLOCKSIZE];          //shared atomic types
  
  int idx = (blockIdx.x*BLOCKSIZE + threadIdx.x);
  short isok = (idx < NAtoms_d);
  idx *= isok;
  
  float  myfor[3], r[3], dist;
  float4 mypos;
  short  mytype;

  mypos  = pos[idx];          //position of the atom for which we will calculate the potential contribution
  mytype = AtomType_d[idx];   //this atom type

  myfor[0] = 0.0f; myfor[1] = 0.0f; myfor[2] = 0.0f;  //set the force to 0

  //printf("the i loop goes to %i\n",NAtoms_d/BLOCKSIZE + 1);

  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL elements!!!)
  for(int i=0; i < NAtoms_d/BLOCKSIZE + 1; i++)
    {

      //every thread loads one position/charge & types from the global memory
      if(threadIdx.x + i*BLOCKSIZE < NAtoms_d)
	{
	  spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE];
	  types[threadIdx.x] = AtomType_d[threadIdx.x + i*BLOCKSIZE];
	}
      __syncthreads();
      //*************************************************************
      
      //now loop on the loaded elements and sum the contributions
      for(int j=0; j<BLOCKSIZE; j++)
	{

	  //stop summing if the index of j is already out
	  if(i*BLOCKSIZE + j >= NAtoms_d)
	    break;
	  
	  //exclude selfterms
	  if( j+i*BLOCKSIZE != idx )
	    {
	      
	      //compute the distance
	      r[0] = spos[j].x - mypos.x;
	      r[1] = spos[j].y - mypos.y;
	      r[2] = spos[j].z - mypos.z;
	      dist = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	      dist = sqrtf(dist);
	      
	      //get the coulomb part
	      //make sure we do not count self shell interactions
	      CoulombForce(r, dist, spos[j].w*mypos.w, myfor );
	      
	      //pair potential part
	      PairForce(r, dist, mytype, types[j], myfor);
	      
	    }
	}

      __syncthreads();
      
    }

  //save the results in the global memory
  if(isok == 1)
    forces[idx] = make_float3(myfor[0],myfor[1],myfor[2]);
  
}















/*
__global__ void Sort_Forces( float *Springs_d, float3 *forces )
{

  int idx = blockIdx.x*BLOCKSIZE + threadIdx.x;
  
  short isok = 1;

  //prevent overflow
  if(idx >= NAtoms_d)
    {
      idx = 0;
      isok = 0;
    }
 
  //if there is no shell...
  if( (Springs_d[idx] == 0.0f) && (isok == 1) )
    {
      forces[idx].x += forces[idx+NAtoms_d].x; //dump the shell force on the core
      forces[idx].y += forces[idx+NAtoms_d].y; //dump the shell force on the core
      forces[idx].z += forces[idx+NAtoms_d].z; //dump the shell force on the core
      forces[idx+NAtoms_d] = make_float3(0.0f,0.0f,0.0f);
    }


}
*/
