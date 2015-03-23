// -*- c++ -*-



//computes the coulomb forces between a and b
__device__ void CoulombForce(float4 a, float4 b, float Force[])
{
  //float3 Force;
  float r[3];
  
  
  r[0] = (b.x - a.x)/0.529177f;
  r[1] = (b.y - a.y)/0.529177f;
  r[2] = (b.z - a.z)/0.529177f;

  float distsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  float dist = rsqrtf(distsq);
  float qq = a.w*b.w*dist / distsq;

  #pragma unroll
  for(int c=0;c<3;c++)
    {
      //r[c] *= dist;
      Force[c] -= r[c]*qq*51.422107f ;
    } 
  
}




__global__ void Calc_Coulomb_Forces( float4 *pos, float3 *forces )
{
  
  __shared__ float4 spos[BLOCKSIZE];
  
  int idx = blockIdx.x*BLOCKSIZE + threadIdx.x;
  float4 mypos;
  float myforce[3];
  short isok = 1;

  if(idx >= 2*NAtoms_d)
    {
      idx = 0;
      isok = 0;
    }
    
  mypos = pos[idx];   //position of the atom for which we will calculate the forces
  myforce[0] = 0.0f;  myforce[1] = 0.0f;  myforce[2] = 0.0f;   //set the forces to 0
  
  
  //loop on all groups of consecutive BLOCKSIZE particles
  for(int i=0; i < 2*NAtoms_d/BLOCKSIZE + 1; i++)
    {
      //every thread loads one position/charge from the global memory
      if(threadIdx.x + i*BLOCKSIZE < 2*NAtoms_d)
	spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE];

      __syncthreads();   //wait for all to finish loading

      //now loop on the loaded elements
      for(int j=0; j<BLOCKSIZE; j++)
	{
 	      
	  //makes sure that we do not have self interaction
	  //nor interaction with itz own shell/core
	  if( (j+i*BLOCKSIZE == idx) || (j+i*BLOCKSIZE == idx+NAtoms_d) ||
	      (j+i*BLOCKSIZE == idx-NAtoms_d) )
	    continue;
	      
	  if(i*BLOCKSIZE + j >= 2*NAtoms_d)
	    break;
	      
	      
	  CoulombForce(mypos, spos[j], myforce);
	  
	}
      __syncthreads();
    }
  
  //copy back from shared to global
  if(isok == 1)
    forces[idx] = make_float3(myforce[0],myforce[1],myforce[2]);

  
}

