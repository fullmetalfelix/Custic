// -*- c++ -*-


//*** ***** FORCES EVALUATION KERNELS (MULTI IMAGE!) ***** ***//
//***********************************************
//no energy evaluation... only forces for many images

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

__global__ void Calc_Both( float4 *pos, short *AtomType_d, float3 *forces )
{
  
  __shared__ float4 spos[BLOCKSIZE];           //shared positions
  __shared__ short  types[BLOCKSIZE];          //shared atomic types
  
  float  myfor[3], r[3], dist; // 7 floats
  float4 mypos;                // 11
  short  mytype;                 
  
  //get the right block index (separation between force and energy)
  int blockx = blockIdx.x;


  //get the image index
  int imgidx = (blockx/BperIMG_d); //index of my image in shared to save space
  blockx -= imgidx*BperIMG_d;

  int idx = (blockx * BLOCKSIZE + threadIdx.x); //index of the atom this thread will take care of
  short isok = (idx < NAtoms_d);
  idx *= isok;

  int windex = idx + imgidx*NAtoms_d;  //16 floats used in registers!!!


  mypos  = pos[windex];       //position of the atom for which we will calculate the potential contribution
  mytype = AtomType_d[idx];//this atom type
  myfor[0] = 0.0f;
  myfor[1] = 0.0f;
  myfor[2] = 0.0f;  //set the force to 0


  //return;

  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL atoms!!!)
  for(int i=0; i < BperIMG_d; i++) //loops over the Blocks needed for one image of the system (NAtoms!)
    {

      //every thread loads one position/charge & types from the global memory
      if(threadIdx.x + i*BLOCKSIZE < NAtoms_d)
	{
	  spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE+NAtoms_d*imgidx];
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
	      CoulombForce(r, dist, spos[j].w*mypos.w, myfor );
	      //pair potential part
	      PairForce(r, dist, mytype, types[j], myfor);
	      
	    }
	}
      //-------------------------------------------------------------


      __syncthreads();
      
    }

  //save the results in the global memory
  //__syncthreads();
  if(isok == 1){
    forces[windex] = make_float3(myfor[0],myfor[1],myfor[2]);
    //forces[windex] = make_float3(blockIdx.x*100.0f+threadIdx.x, blockx,windex);
  } 

}



int GetForcesEnergy()
{
  dim3 dimBlock_f(BLOCKSIZE);
  dim3 dimGrid_f(Run.BperIMG * Run.MDimages);
  if(dimGrid_f.x == 0) dimGrid_f = 1;
  
  //int i;
  unsigned int t1;

  cutCreateTimer(&t1);
  cutStartTimer(t1);

  Calc_Both<<<dimGrid_f, dimBlock_f>>>(Charges_d, AtomType_d, forces_d);
  cudaThreadSynchronize();
  
  cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);           //copy back the forces
  
  cutStopTimer(t1);

  printf("forces time: %f ms (witho copyback)\n",cutGetTimerValue(t1));
  cutDeleteTimer(t1);
  
    FILE *fp = fopen("forceold.out","w");
    
    //for(img=0;img<Run.MDimages;img++)
    for(int i=0;i<Run.MDimages*Run.NAtoms;i++)
    {
    //fprintf(fp,"%8.5f %8.5f %8.5f - %8.5f \n",Charges_h[i].x,Charges_h[i].y,Charges_h[i].z,Charges_h[i].w);
    fprintf(fp,"%8.5f %8.5f %8.5f \n",forces_h[i].x,forces_h[i].y,forces_h[i].z);
    //fprintf(fp,"atom%i %f %f %f ...  %f\n",i,forces_h[i].x-forces_h[i+Run.NAtoms].x, forces_h[i].y-forces_h[i+Run.NAtoms].y, 
    //	forces_h[i].z-forces_h[i+Run.NAtoms].z,Charges_h[i].w);
    //fprintf(fp,"index%i blk%f  thr%f  aidx%f  isok %f\n",i,forces_h[i].x, forces_h[i].y,forces_h[i].z,
    //	  Potentials_h[i]);
    //fprintf(fp,"index%i %f  %f %f %f\n",i,Potentials_h[i],Charges_h[i].x,Charges_h[i].y,Charges_h[i].z);
    }
    fclose(fp);
  
  
  return true;
  
}
