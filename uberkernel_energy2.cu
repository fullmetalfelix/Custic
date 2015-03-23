// -*- c++ -*-

extern "C"
{
//*** ***** ENERGY EVALUATION KERNELS (MULTI IMAGE) ***** ***//
//***********************************************


__device__ float PairPotential(float dist, int i, int j)
{
  float4 p =  ShellParameters_d[ max(i,j)*(max(i,j)+1)/2 + min(i,j) ];  //!! CONST MEMORY READ

  return p.x*expf(-dist/p.y) + p.z*powf(dist,p.w);
}


__global__ void Calc_Energy( float4 *pos, short *AtomType_d, float *Potentials_d )
{
  
  __shared__ float4 spos[BLOCKSIZE];           //shared positions
  __shared__ short  types[BLOCKSIZE];          //shared atomic types
  
  float  r[3], dist;           // 4 floats
  float  mypot = 0.0f;         // 5
  float4 mypos;                // 9
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


  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL atoms!!!)
  for(int i=0; i < BperIMG_d; i++){ //loops over the Blocks needed for one image of the system (NAtoms!)
    
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
	      
	      
	    if(isE)
	      {
		mypot +=  0.5f * spos[j].w*mypos.w * 14.39964524f / dist;
		mypot +=  0.5f * PairPotential( dist, mytype, types[j] );	    
	      }
	    else
	      {
		//get the coulomb part
		CoulombForce(r, dist, spos[j].w*mypos.w, myfor );
		//pair potential part
		PairForce(r, dist, mytype, types[j], myfor);
	      }
	      
	  }
      }
    //-------------------------------------------------------------
    __syncthreads();
      
  }

  //save the results in the global memory
  spos[threadIdx.x].x = mypot * isok; //every thread writes the partial in the shared
  __syncthreads();

  //the first thread of each block copies the result in the global memory
  mypot = 0.0f;
  if(threadIdx.x == 0)
    {
      for(int i=0;i<BLOCKSIZE;i++)
	  mypot += spos[i].x;
      
      Potentials_d[blockx+BperIMG_d*imgidx] = mypot;
    }
  
}


//get the energy of the system
int GetForcesEnergy()
  {
    dim3 dimBlock_f(BLOCKSIZE);
    dim3 dimGrid_f(Run.BperIMG * Run.MDimages);
    if(dimGrid_f.x == 0) dimGrid_f = 2;
    
    int i,img;
    float energy;
    
    
    Calc_Energy<<<dimGrid_f, dimBlock_f>>>(Charges_d, AtomType_d, Potentials_d);
    cudaThreadSynchronize();
    
    cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);           //copy back the forces
    cudaMemcpy(Potentials_h, Potentials_d, f1_NObj, cudaMemcpyDeviceToHost);

    //sum up the potentials given by each block
    for(img=0; img<Run.MDimages; img++)
      {
	energy = 0.0f;
	for(i=0;i<Run.BperIMG;i++)
	  energy += Potentials_h[i+img*Run.BperIMG];
	Etot[img] = energy;
	//printf("Energy[%i] is %f\n",img,energy);
      }
   


    
    return true;

  }


}
