// -*- c++ -*-

//MD functions
extern "C"
{

  // Kb = 1.3806503e-23 SI
  // 1 uma = 1.660538782e-27 Kg
  // 1 eV  = 1.60217646e-19 J
  // 1 eV/(Åuma) = 0.964853382e18  m/s² = 0.964853382e-2 Å/fs² (acceleration)
  // 1 N = 0.602214179e7 uma Å/fs²
  // 1 Å/fs = 100000 m/s

  //positions are in Å
  //speeds are in Å/fs
  //forces are in eV/Å
  //time is in fs
  //masses are in uma


  __device__ float3 operator+(const float3 &a, const float3 &b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
  }
  __device__ float3 operator*(const float &a, const float3 &b) {
    return make_float3(a*b.x, a*b.y, a*b.z);
  }

  //this is used to do pos-holderpos and put the spring (holder.w) in the fourth component
  __device__ float4 operator-(const float4 &pos, const float4 &holder) {
    return make_float4(pos.x-holder.x, pos.y-holder.y, pos.z-holder.z,holder.w);
  }

  

  /* Sum forces on the tip:
     do a parallel reduction to find forces on the tip; all images, one big block per image.
     MaxTip: number of tip atoms, it assumes the tip atoms are aligned in the beginning of the list,
     if not, you will die!
     
  *//*
  __global__ void GetTipForce_GPU_v2(float3 *Force, float3 *TipF, float4 *Pos, float4 *Holder)
  {
    __shared__ float3 sForce[BLOCKFAST];

    float4 sHolder;
    
    float3 tmpF,finF;
    short offset,j,k;
    //short img = blockIdx.x*NAtoms_d; //one block per image
    
    finF = TipF[blockIdx.x];
    
    for(j=0; j<ceilf((float)NTipAtoms_d/BLOCKFAST); j++){
      
      //precache with offset
      sForce[threadIdx.x] = make_float3(0.0f,0.0f,0.0f);
      if(threadIdx.x + BLOCKFAST*j < NTipAtoms_d){
	sForce[threadIdx.x] = Force[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j];
	sHolder = Pos[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j]-
	          Holder[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j];

	if(sHolder.w > 0){
	  sForce[threadIdx.x].x += sHolder.w * sHolder.x; //remove the spring force
	  sForce[threadIdx.x].y += sHolder.w * sHolder.y;
	  sForce[threadIdx.x].z += sHolder.w * sHolder.z;
	}
      }
      __syncthreads();
      //--------------------
      
      //sum the fetched guys
      for(offset = 1; offset < BLOCKFAST; offset *= 2){
	
	k = (threadIdx.x+offset);  //j is the second addend
	tmpF = make_float3(0.0f,0.0f,0.0f); //reset the temp
	
	if ( k < BLOCKFAST ) //if not going out of bounds
	  tmpF = sForce[k];  //set it to the next addend
	
	tmpF = sForce[threadIdx.x] + tmpF; //so this adds 0 to threads targeting OOB
	__syncthreads();
      
	sForce[threadIdx.x] = tmpF; //swap the registers
	__syncthreads();
      }
      
      if(threadIdx.x == 0){ //thread 0 stores the result
	finF = finF + tmpF;
      }
      
    }

    //now thread 0 has the final force
    if(threadIdx.x == 0)
      TipF[blockIdx.x] = finF;

  }
    */
  //WHAT if forces would be float4 instead?? we'll never know!


  //ONLY NORMAL FORCE!
  __global__ void GetTipForce_GPU_v3(float3 *Force, float3 *TipF, float4 *Pos, float4 *Holder)
  {
    __shared__ float sForce[BLOCKFAST];

    float4 sHolder;
    
    float tmpF,finF;
    short offset,j,k;
    //short img = blockIdx.x*NAtoms_d; //one block per image
    
    finF = TipF[blockIdx.x].z;
    
    for(j=0; j<ceilf((float)NTipAtoms_d/BLOCKFAST); j++){
      
      //precache with offset
      sForce[threadIdx.x] = 0.0f;
      if(threadIdx.x + BLOCKFAST*j < NTipAtoms_d){
	sForce[threadIdx.x] = Force[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j].z;
	sHolder = Pos[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j]-
	          Holder[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j];

	if(sHolder.w > 0){
	  sForce[threadIdx.x] += sHolder.w * sHolder.z;//remove the spring force
	}
      }
      __syncthreads();
      //--------------------
      
      //sum the fetched guys
      for(offset = 1; offset < BLOCKFAST; offset *= 2){
	
	k = (threadIdx.x+offset);  //j is the second addend
	tmpF = 0.0f; //reset the temp
	
	if ( k < BLOCKFAST ) //if not going out of bounds
	  tmpF = sForce[k];  //set it to the next addend
	
	tmpF = sForce[threadIdx.x] + tmpF; //so this adds 0 to threads targeting OOB
	__syncthreads();
      
	sForce[threadIdx.x] = tmpF; //swap the registers
	__syncthreads();
      }
      
      if(threadIdx.x == 0){ //thread 0 stores the result
	finF += tmpF;
      }
      
    }

    //now thread 0 has the final force
    if(threadIdx.x == 0)
      TipF[blockIdx.x].z = finF;

  }

  //Lateral ONLY????????????????????????????++NO! this does all components!
  __global__ void GetTipForce_GPU_vlat(float3 *Force, float3 *TipF, float4 *Pos, float4 *Holder)
  {
    __shared__ float3 sForce[BLOCKFAST];

    float4 sHolder;
    
    float3 tmpF,finF;
    short offset,j,k;
    //short img = blockIdx.x*NAtoms_d; //one block per image
    
    finF = TipF[blockIdx.x];
    
    for(j=0; j<ceilf((float)NTipAtoms_d/BLOCKFAST); j++){
      
      //precache with offset
      sForce[threadIdx.x] = make_float3(0.0f, 0.0f, 0.0f);
      if(threadIdx.x + BLOCKFAST*j < NTipAtoms_d){
	sForce[threadIdx.x] = Force[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j];
	sHolder = Pos[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j]-
	          Holder[blockIdx.x*NAtoms_d+threadIdx.x + BLOCKFAST*j];

	if(sHolder.w > 0){
	  sForce[threadIdx.x].x += sHolder.w * sHolder.x;//remove the spring force
	  sForce[threadIdx.x].y += sHolder.w * sHolder.y;//remove the spring force
	  sForce[threadIdx.x].z += sHolder.w * sHolder.z;//remove the spring force
	}
      }
      __syncthreads();
      //--------------------
      
      //sum the fetched guys
      for(offset = 1; offset < BLOCKFAST; offset *= 2){
	
	k = (threadIdx.x+offset);  //j is the second addend
	tmpF = make_float3(0.0f,0.0f,0.0f); //reset the temp
	
	if ( k < BLOCKFAST ) //if not going out of bounds
	  tmpF = sForce[k];  //set it to the next addend
	
	tmpF = sForce[threadIdx.x] + tmpF; //so this adds 0 to threads targeting OOB
	__syncthreads();
      
	sForce[threadIdx.x] = tmpF; //swap the registers
	__syncthreads();
      }
      
      if(threadIdx.x == 0){ //thread 0 stores the result
	finF = finF + tmpF;
      }
      
    }

    //now thread 0 has the final force
    if(threadIdx.x == 0)
      TipF[blockIdx.x] = finF;

  }
  
  void GetTipForce(  )
  {
    GetTipForce_GPU_vlat<<<Run.MDimages, FastBlock>>>(forces_d, TipForce_GPU, Charges_d, Holders_d);
    cudaThreadSynchronize();
  }
  

  //DEBUG VERSION!
  /*
  void GetTipForce_cpu( void );

  void GetTipForce(  )
  {
    //test routine
        
    unsigned int t1,t2;
    
   
    cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);
   
    //time the CPU calculation
    cutCreateTimer(&t1);
    cutStartTimer(t1);
    
    for(int i=0;i<100;i++){
      cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost); //time the transfer too
      GetTipForce_cpu( );
    } 

    cutStopTimer(t1);
    //----------------------------------

    
    //---GPU calculation---
    dim3 grid(Run.MDimages); //one block per image
    dim3 blok(BLOCKFAST);

    float3 *gpuf;

    gpuf = (float3*)malloc(Run.MDimages*sizeof(float3));
    for(int i=0;i<Run.MDimages;i++)
      gpuf[i] = make_float3(0.0f,0.0f,0.0f);
    cudaMemcpy(TipForce_GPU, gpuf, Run.MDimages*sizeof(float3), cudaMemcpyHostToDevice);

    cutCreateTimer(&t2);cutStartTimer(t2);
    
    for(int i=0;i<100;i++){
      GetTipForce_GPU_v1<<<FastGrid, blok>>>(forces_d, TipAtoms_d, TipForce_GPU, BlockForce_GPU, ForceDone_GPU);
      cudaThreadSynchronize();
    }
  
    cutStopTimer(t2);
 
    cudaMemcpy(gpuf, TipForce_GPU, Run.MDimages*sizeof(float3),cudaMemcpyDeviceToHost );
    cudaThreadSynchronize();
    
    printf("time: %f - %f ms\n",cutGetTimerValue(t1)/100.0f, cutGetTimerValue(t2)/100.0f);

    for(int i=0;i<Run.MDimages;i++)
      printf("result: %f - %f \n",ClusterTemps[i],gpuf[i].z);  
        

    free(gpuf);
    cutDeleteTimer(t1);
    cutDeleteTimer(t2);
  }




 //calculate the temperatures of the clusters - optimized
  void GetTipForce_cpu(  )
  {
    int i,j;
    float hv;
    float3 F[Run.MDimages];    //cluster temp
    float3 frc;

    for(j=0; j<Run.MDimages; j++){
      F[j] = make_float3(0.0f,0.0f,0.0f);
      
      for(i=0;i<Run.NAtoms;i++){  //sum all kins
	if(TipAtoms_h[i] == 1){
	  frc = forces_h[i+j*Run.NAtoms];
	  F[j].x += frc.x;  
	  F[j].y += frc.y;  
	  F[j].z += frc.z;
	}
      }
      ClusterTemps[j] += F[j].z;
    }
  }
  */

}


 
