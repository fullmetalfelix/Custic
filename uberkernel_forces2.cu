// -*- c++ -*-


extern "C"
{
//*** ***** FORCES EVALUATION KERNELS (MULTI IMAGE!) ***** ***//
//***********************************************
//no energy evaluation... only forces for many images

  __device__ float4 operator-(const float4 &pos, const float4 &holder);
__global__ void Calc_Both_v7p( float4 *pos, short *AtomType_d, float3 *forces, float4 *hld )
{
  
  __shared__ float4 spos[BLOCKSIZE];           //shared positions
  __shared__ short  types[BLOCKSIZE];          //shared atomic types
  __shared__ float3 spots[MAX_PPOT];


  float3 myfor;
  float4 r;
  float4 mypos;                // 11              
  float3 p;


  short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
  asd.x = (blockIdx.x/BperIMG_d);
  asd.y = blockIdx.x - asd.x*BperIMG_d;
  asd.z = asd.y * BLOCKSIZE + threadIdx.x;
  asd.x *= NAtoms_d;  // x=imageindex*Natoms

  short isok = (asd.z < NAtoms_d);
  asd.z *= isok;


  mypos = pos[asd.z+asd.x];       //position of the atom for which we will calculate the potential contribution
  asd.w = AtomType_d[asd.z];//this atom type
  myfor.x = 0.0f;
  myfor.y = 0.0f;
  myfor.z = 0.0f;  //set the force to 0
  
  float qq;
  
  //prefetch pair pots
  if(threadIdx.x < MAX_PPOT){
    spots[threadIdx.x] = ShellParameters_d[threadIdx.x];
  }
  __syncthreads();
  

  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL atoms!!!)
  for(int i=0; i < BperIMG_d; i++) //loops over the Blocks needed for one image of the system (NAtoms!)
    {

      //every thread loads one position/charge & types from the global memory
      if(threadIdx.x + i*BLOCKSIZE < NAtoms_d)
	{
	  spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE+asd.x];
	  types[threadIdx.x] = AtomType_d[threadIdx.x + i*BLOCKSIZE];
	}
      __syncthreads();
      //*************************************************************
      
      //now loop on the loaded elements and sum the contributions
      for(int j=0; j<BLOCKSIZE; j++){
	
	//stop summing if the index of j is already out
	if(i*BLOCKSIZE + j >= NAtoms_d)
	  break;
	
	//exclude selfterms
	if( j+i*BLOCKSIZE != asd.z ){
	  
	  //compute the distance
	  r.x = spos[j].x - mypos.x;
	  r.y = spos[j].y - mypos.y;
	  r.z = spos[j].z - mypos.z;
	  r.w = r.x*r.x + r.y*r.y + r.z*r.z; //r.w = r**2
	  p.x = r.w;
	  r.w = sqrtf(r.w);   //r.w = |r|
	  
	  //get the coulomb part
	  p.x = r.w * p.x; //p.x = r**3
	    //qq = r.w * p.x; //<-- r**3
	  p.y = spos[j].w*mypos.w *14.3996368420f / p.x;
	  
	  myfor.x -= r.x*p.y;
	  myfor.y -= r.y*p.y;
	  myfor.z -= r.z*p.y;

	  //pair pot
	  qq = -6.0f/(p.x*p.x*r.w); // = r**-7

	  p = spots[asd.w*MAXTYPES+types[j]];
	  qq = (-p.x*expf(-r.w*p.y)*p.y + p.z*qq) / r.w;
	  
	  //if(threadIdx.x == 0)
	  //printf("pair interacton %f\n",qq*r[0]);
	  
	  myfor.x += r.x*qq;
	  myfor.y += r.y*qq;
	  myfor.z += r.z*qq;
	}
      }
      //-------------------------------------------------------------
      
      
      __syncthreads();
      
    }
  
  //add the spring part
  //spos[threadIdx.x] = hld[asd.z]; //the holder is read from the first image!!! (they are all the same!)
  
  mypos = mypos - hld[asd.z];
  if( mypos.w > 0 ){
    myfor.x -= mypos.w * mypos.x;
    myfor.y -= mypos.w * mypos.y;
    myfor.z -= mypos.w * mypos.z;
  }
  //force is 0 on the fixed atoms
  /*
  if( spos[threadIdx.x].w < 0.0f ){
    myfor.x = 0.0f;
    myfor.y = 0.0f;
    myfor.z = 0.0f;
    }*/

  //save the results in the global memory
  if(isok == 1)
    forces[asd.x+asd.z] = myfor;

}


int GetForcesEnergy()
{
  
  Calc_Both_v7p<<<SlowGrid, SlowBlock>>>(Charges_d, AtomType_d, forces_d, Holders_d);
  cudaThreadSynchronize();
  return true;
  
}


 /*
// DEBUG VERSION!
int GetForcesEnergy()
{
  //dim3 dimBlock_f(BLOCKSIZE);
  //dim3 dimGrid_f(Run.BperIMG * Run.MDimages);
  //if(dimGrid_f.x == 0) dimGrid_f = 1;
  unsigned int t1;

  cutCreateTimer(&t1);
  cutStartTimer(t1);

  //int i;
  for(int i=0;i<1000;i++){
    Calc_Both_v7p<<<SlowGrid, SlowBlock>>>(Charges_d, AtomType_d, forces_d, Holders_d);
    cudaThreadSynchronize();
  }
   
  cutStopTimer(t1);

  printf("forces time: %f ms (witho copyback)\n",cutGetTimerValue(t1)/1000);
  cutDeleteTimer(t1);

  cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);           //copy back the forces
  
    FILE *fp = fopen("force1.out","w");
    
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
 */

}
