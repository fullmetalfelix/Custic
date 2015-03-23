// -*- c++ -*-


//*** ***** FORCES EVALUATION KERNELS (MULTI IMAGE!) ***** ***//
//***********************************************
//no energy evaluation... only forces for many images

__device__ float3 CoulombForce(float4 r, float Cc, float3 Force)
{

  //dist is in r.w
  float qq = r.w * r.w * r.w;
  qq = Cc *14.3996368420f / qq;

  
  Force.x -= r.x*qq;
  Force.y -= r.y*qq;
  Force.z -= r.z*qq;
  return Force;
}

__device__ float3 PairForce(float4 r, int i, int j, float3 Force)
{

  float4 p =  ShellParameters_d[ max(i,j)*(max(i,j)+1)/2 + min(i,j) ];  //!! CONST MEMORY READ

  //this is -grad(E)
  float qq = (-p.x*expf(-r.w/p.y)/p.y + p.z*p.w*powf(r.w,p.w-1)) / r.w;
 
  Force.x += r.x*qq;
  Force.y += r.y*qq;
  Force.z += r.z*qq;

  return Force;
 
}

__global__ void Calc_Coulomb( float4 *pos, float3 *forces )
{
  
  __shared__ float4 spos[BLOCKSIZE];           //shared positions  (64threads => 1024 bytes)

  //__shared__ float4 mypos[BLOCKSIZE];
  //__shared__ float4 myhld[BLOCKSIZE];
  //__shared__ float3 myfor[BLOCKSIZE];
  
  float4 mypos, r;
  float3 myfor;

  //float  dist; // 7 floats
  //float4  myhld;         // 15
  //short  mytype;                 
  
  //get the right block index (separation between force and energy)
  short3 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
  asd.x = (blockIdx.x/BperIMG_d);
  asd.y = blockIdx.x - asd.x*BperIMG_d;
  asd.z = asd.y * BLOCKSIZE + threadIdx.x;
  //asd.w = AtomType_d[asd.z];
  
  float qq;

  short isok = (asd.z < NAtoms_d);
  asd.z *= isok;

  //int windex = idx + imgidx*NAtoms_d;  //16 floats used in registers!!!


  mypos  = pos[asd.z + NAtoms_d*asd.x];       //position of the atom for which we will calculate the potential contribution
  //mytype = AtomType_d[asd.z];//this atom type
  myfor.x = 0.0f;
  myfor.y = 0.0f;
  myfor.z = 0.0f;  //set the force to 0


  //return;

  //loop on all groups of consecutive BLOCKSIZE particles (will loop AAALL atoms!!!)
  for(int i=0; i < BperIMG_d; i++) //loops over the Blocks needed for one image of the system (NAtoms!)
    {

      //every thread loads one position/charge & types from the global memory
      if(threadIdx.x + i*BLOCKSIZE < NAtoms_d){
	spos[threadIdx.x] = pos[threadIdx.x + i*BLOCKSIZE + NAtoms_d*asd.x];
	//types[threadIdx.x] = AtomType_d[threadIdx.x + i*BLOCKSIZE];
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
	  if( j+i*BLOCKSIZE != asd.z )
	    {
	      
	      //compute the distance
	      r.x = spos[j].x - mypos.x;
	      r.y = spos[j].y - mypos.y;
	      r.z = spos[j].z - mypos.z;
	      r.w = r.x*r.x + r.y*r.y + r.z*r.z;
	      r.w = sqrtf(r.w);
	      
	      //get the coulomb part
	      //myfor = CoulombForce(r, spos[j].w*mypos.w, myfor );
	      qq = r.w * r.w * r.w;
	      qq = spos[j].w*mypos.w *14.3996368420f / qq;

	      myfor.x -= r.x*qq;
	      myfor.y -= r.y*qq;
	      myfor.z -= r.z*qq;



	      //pair potential part
	      //myfor = PairForce(r, asd.w, types[j], myfor);
	      /*
	      p =  ShellParameters_d[ max(asd.w,types[j])*(max(asd.w,types[j])+1)/2 + min(asd.w,types[j]) ];  
	      //!! CONST MEMORY READ

	      //this is -grad(E)
	      qq = (-p.x*expf(-r.w/p.y)/p.y + p.z*p.w*powf(r.w,p.w-1)) / r.w;
 
	      myfor.x += r.x*qq;
	      myfor.y += r.y*qq;
	      myfor.z += r.z*qq;*/


	      
	    }
	}
      //-------------------------------------------------------------


      __syncthreads();
      
    }

  //save the results in the global memory
  //__syncthreads();
  if(isok == 1){
    forces[asd.z + NAtoms_d*asd.x] = myfor;
  } 

}



int GetForcesEnergy()
{
  dim3 dimBlock_f(BLOCKSIZE);
  dim3 dimGrid_f(Run.BperIMG * Run.MDimages);
  if(dimGrid_f.x == 0) dimGrid_f = 1;
  unsigned int t1;

  cutCreateTimer(&t1);
  cutStartTimer(t1);

  //int i;
  
  Calc_Coulomb<<<dimGrid_f, dimBlock_f>>>(Charges_d, forces_d);
  cudaThreadSynchronize();
  
  cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);           //copy back the forces
  
  cutStopTimer(t1);

  printf("forces time: %f ms (witho copyback)\n",cutGetTimerValue(t1));
  cutDeleteTimer(t1);

  
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
