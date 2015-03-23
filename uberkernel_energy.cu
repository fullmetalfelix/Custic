// -*- c++ -*-

extern "C"
{
//*** ***** ENERGY EVALUATION KERNELS ***** ***//
//***********************************************


__device__ float PairPotential(float dist, int i, int j)
{
  float4 p =  ShellParameters_d[ max(i,j)*(max(i,j)+1)/2 + min(i,j) ];  //!! CONST MEMORY READ

  return p.x*expf(-dist/p.y) + p.z*powf(dist,p.w);
}


__global__ void Calc_Energy( float4 *pos, short *AtomType_d, float *Potentials_d )
{
  
  __shared__ float4 spos[BLOCKSIZE];
  __shared__ float  spots[BLOCKSIZE];
  __shared__ short  types[BLOCKSIZE];


  int idx = blockIdx.x*BLOCKSIZE + threadIdx.x;
  short isok = (idx < NAtoms_d);
  idx *= isok;

  int i;
  float  mypot = 0.0f;
  float r[3], dist;
  float4 mypos;
  short  mytype;

 
  mypos  = pos[idx];   //position of the atom for which we will calculate the potential contribution
  mytype = AtomType_d[idx]; 

  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL elements!!!)
  for(i=0; i < NAtoms_d/BLOCKSIZE + 1; i++)
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
	      mypot +=  0.5f * spos[j].w*mypos.w * 14.39964524f / dist;
	      //if(threadIdx.x == 0)
	      //printf("ADD Coulomb %i %i: %f\n",i,j,0.5f * spos[j].w*mypos.w * 14.39964524f / dist);
	      
	      //pair potential part
	      //do this only between shells take only shells
	      mypot +=  0.5f * PairPotential( dist, mytype, types[j] );	      

	    }
	}

      __syncthreads();
      
    }

   

  //save the results in the shared memory
  spots[threadIdx.x] = mypot * isok;
  __syncthreads();

  //the first thread of each block copies the result in the global memory
  mypot = 0.0f;
  if(threadIdx.x == 0)
    {
      for(i=0;i<BLOCKSIZE;i++)
	{
	  mypot += spots[i];
	  //printf("spots[%i] = %f\n",i,spots[i]);
	}
      Potentials_d[blockIdx.x] = mypot;
    }
  
  
}


//get the energy of the system
float GetEnergy(int DebugPrint)
{
  
  //unsigned int timer1;
  
  dim3 dimBlock_f(BLOCKSIZE);
  dim3 dimGrid_f(Run.NAtoms / BLOCKSIZE + 1);
  
  int i;
  float energy = 0.0f;
  
  for(i=0;i<Run.NAtoms;i++)
    Potentials_h[i] = 0.0f;
  
  //do the kernel
  Calc_Energy<<<dimGrid_f, dimBlock_f>>>(Charges_d, AtomType_d, Potentials_d);
  cudaThreadSynchronize();
  
  cudaMemcpy(Potentials_h, Potentials_d, f1_NObj, cudaMemcpyDeviceToHost);
  
  //sum up the potentials given by each block
  for(i=0;i<dimGrid_f.x;i++)
    energy += Potentials_h[i];
  
  //printf("ENERGY!\n");
  
  return energy;
}



}
