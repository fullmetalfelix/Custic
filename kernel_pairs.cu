// -*- c++ -*-



//ONLY FOR BUCKINGHAM TYPE!
__device__ void PairForce(float4 a, float4 b, float4 p, float Force[])
{

  float r[3];
  
  
  r[0] = b.x - a.x;
  r[1] = b.y - a.y;
  r[2] = b.z - a.z;

  float distsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  float rdist = rsqrtf(distsq);
  float dist = 1.0f/rdist; //sqrtf(distsq);

  
  float qq = +p.x*expf(-dist/p.y)/p.y + p.z*p.w*powf(dist,p.w-1);

  #pragma unroll
  for(int c=0;c<3;c++)
    {
      r[c] *= rdist;
      Force[c] += r[c]*qq ;
    } 


}


__global__ void Calc_Pairs( float4 *pos, float3 *forces, float *Springs_d, short *AtomType_d, float4 *ShellParameters_d)
{
  
  __shared__ float4 spos[BLOCKSIZE];   //shared positions
  __shared__ short stypes[BLOCKSIZE];    //shared atom types
  __shared__ float4 shellpars[64];     //shell-shell parameters (1 KiB)

  int i,j;
  int idx = blockIdx.x*BLOCKSIZE + threadIdx.x + NAtoms_d; //start from the first shell!
  short isok = 1;

  float myforce[3]; //resulting force
  float4 mypos;
  short myType;

  //load all shell parameters in the shared memory
  for(i=0; i < NTypes_d*NTypes_d/BLOCKSIZE + 1; i++)
  {
      if(threadIdx.x + i*BLOCKSIZE < NTypes_d*NTypes_d)   //if the thread will get a valid element...
          shellpars[threadIdx.x + i*BLOCKSIZE] = ShellParameters_d[threadIdx.x + i*BLOCKSIZE];
  }
  __syncthreads();


  //reset the index for outofbounds shit
  if(idx >= 2*NAtoms_d)
    {
      idx = NAtoms_d;
      isok = 0;
    }

  mypos = pos[idx];
  myType = AtomType_d[idx-NAtoms_d];


  //set the object index: if there is no shell, act on the core
  int ObjIndex;
  if(Springs_d[idx-NAtoms_d] == 0.0f)
    ObjIndex = idx - NAtoms_d; //use the index of the core
  else
    ObjIndex = idx;            //or the index of the shell.

  //set the force to 0
  myforce[0]=0.0f; myforce[1]=0.0f; myforce[2]=0.0f;
  

  //loop on all groups of consecutive BLOCKSIZE particles
  for(i=0; i<NAtoms_d/BLOCKSIZE + 1; i++)
    {
      
      //every thread loads one position/charge from the global memory (no OOB allowed)
      if(threadIdx.x + i*BLOCKSIZE <  NAtoms_d)
	{
	  spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE + NAtoms_d]; //pos,charge
	  stypes[threadIdx.x] = AtomType_d[threadIdx.x + i*BLOCKSIZE];   //type 
	}
      __syncthreads();   //wait for all to finish loading

      //now loop on the loaded elements
      //#pragma unroll
      for(j=0; j<BLOCKSIZE; j++)
	{
          //avoid the self interaction
	  if(idx-NAtoms_d == i*BLOCKSIZE + j)
	    continue;

	  //if we are out of bound stop summing!
	  if(i*BLOCKSIZE + j >= NAtoms_d)
	    break;


	  //compute the short range term
	  PairForce(mypos, spos[j], shellpars[myType*NTypes_d + stypes[j]], myforce);
	  
	  
	}
      __syncthreads();
    }
  
  
  //copy back from shared to global
  //printf("thread %i %i, %f %f %f\n",threadIdx.x,blockIdx.x, myforce[0],myforce[1],myforce[2]);
  if(isok == 1)
    {
      forces[ObjIndex].x -= myforce[0];
      forces[ObjIndex].y -= myforce[1];
      forces[ObjIndex].z -= myforce[2];
    }

  
}
