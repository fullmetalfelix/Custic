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



  /* Calculate the temperature of the WHOLE SYSTEM
     mapping: one thread per atom
     The blocks do partial sums and then use atomics to fetch the other blocks value,
     then thread 0, block 0(referred to the image) finishes the sum and writes back.
     parallel reduction algorithm to sum kinetic energies

     NOTES: this kernel does not take into account fixed atoms!
  */
  __global__ void GetT(float4 *V, float *T)
  {
    __shared__ float TTemp[BLOCKFAST];
    
    float STemp = 0.0f;
    int offset,j;

    float4 spd;

    //set the limits
    short2 lims;
    lims.x = ThermoBounds_d.x;
    lims.y = ThermoBounds_d.y;
    if(blockIdx.y == 1){ //if the block is in second row, calculate the surface
      lims.x = ThermoBounds_d.z;
      lims.y = ThermoBounds_d.w;
    }

    //TTemp[threadIdx.x] = 0.0f; // set shared temp to 0

    for(j=0; j<ceilf((float(lims.y-lims.x))/BLOCKFAST); j++){ //loop over the right amount of blocks
      
      //each thread reads a speed
      //     starts from the limit  skip j blocks  skip blockidx.x images
      offset = threadIdx.x + lims.x + BLOCKFAST*j;
      spd = V[offset + blockIdx.x*NAtoms_d];
      spd.x = spd.x*spd.x + spd.y*spd.y + spd.z*spd.z;
      spd.x *= 51.821345f * spd.w; //kinetic energy in eV!!!
      if( (offset) < lims.y )
	STemp += spd.x;  
    }
    //now all the threads have the kinetic energy of a subset of atoms... do parallel reduction
    TTemp[threadIdx.x] = STemp;
    __syncthreads();

    for(offset = 1; offset < BLOCKFAST; offset *= 2){
      
      j = (threadIdx.x+offset);  //j is the second addend
      STemp = 0.0f; //reset the temp
      
      if ( j < BLOCKFAST ) //if not going out of bounds
	STemp = TTemp[j];  //set it to the next addend
      
      STemp = TTemp[threadIdx.x] + STemp; //so this adds 0 to threads targeting OOB
      __syncthreads();
      
      TTemp[threadIdx.x] = STemp; //swap the registers
      __syncthreads();
    }
    
    if(threadIdx.x == 0){ //       7736.336672065333f
      T[blockIdx.x*2+blockIdx.y] = 7736.336672f * STemp/(lims.y-lims.x);
    }
    
    
  }
  

  
  void GetSystemTemp(  )
  {
    
    GetT<<<TGrid, FastBlock>>>(Speeds_d, Temperature_GPU);
    cudaThreadSynchronize();
    //now each image knows its global temperature
  }
  

  //DEBUG VERSION!
  /*
  void GetClustersTemp_cpu( );
  float *ClusterTemps;
  void GetSystemTemp(  )
  {
    //test routine
    
    dim3 blok(BLOCKFAST);
    dim3 grid(Run.MDimages,2);
    if(grid.x == 0) grid.x = 1;
  

    unsigned int t1,t2;
    
    ClusterTemps = (float*)calloc(2*Run.MDimages,sizeof(float));

    //time the CPU calculation
    cutCreateTimer(&t1);
    cutStartTimer(t1);
    
    for(int i=0;i<100;i++)
      GetClustersTemp_cpu( );

    cutStopTimer(t1);
   
    
    float *gput;
    gput  = (float*)calloc(2*Run.MDimages, sizeof(float));
    //cudaMemcpy(Temperature_GPU,gput,Run.MDimages*sizeof(float),cudaMemcpyHostToDevice );


    cutCreateTimer(&t2);cutStartTimer(t2);
    
    
    for(int i=0;i<100;i++){
      //cudaMemcpy(TempDone_GPU,tdone, Run.MDimages*Run.BperIMG*sizeof(unsigned int), cudaMemcpyHostToDevice);
      GetT<<<grid, blok>>>(Speeds_d, Temperature_GPU);
      cudaThreadSynchronize();
    }
    
    cutStopTimer(t2);
   
    cudaMemcpy(gput,Temperature_GPU, 2*Run.MDimages*sizeof(float), cudaMemcpyDeviceToHost );
    
    printf("time: %f - %f ms\n",cutGetTimerValue(t1)/100.0f, cutGetTimerValue(t2)/100.0f);

    for(int i=0;i<Run.MDimages;i++){
      printf("result: (%f %f) - (%f %f) \n",ClusterTemps[i*2],ClusterTemps[i*2+1],
	     gput[i*2],gput[i*2+1]);  
    }
    
    free(gput); free(ClusterTemps);
    cutDeleteTimer(t1);
    cutDeleteTimer(t2);

  }




 //calculate the temperatures of the clusters - optimized
  void GetClustersTemp_cpu( )
  {
    int i,j;
    float hv;
    float T[2*Run.MDimages];    //cluster temp
    float4 spd;
 
    cudaMemcpy(Speeds_h, Speeds_d, f4_NObj, cudaMemcpyDeviceToHost);

    for(j=0; j<Run.MDimages; j++){
      T[j*2]   = 0.0f;
      T[j*2+1] = 0.0f;
      
      for(i=ThermoBounds_h.x;i<ThermoBounds_h.y;i++){  //sum all kins
	spd = Speeds_h[i+j*Run.NAtoms];
	//printf("%f %f %f\n",spd.x,spd.y,spd.z);
	hv = spd.x*spd.x + spd.y*spd.y + spd.z*spd.z;
	hv = hv * 1.0364269e2 * 0.5f * spd.w; //kinetic energy in eV!!!
	T[j*2] += hv;  
      }
      ClusterTemps[j*2] = 2.0f*T[j*2]/(3.0f*(ThermoBounds_h.y-ThermoBounds_h.x)) * 11604.505008098f;
      
      for(i=ThermoBounds_h.z;i<ThermoBounds_h.w;i++){  //sum all kins
	spd = Speeds_h[i+j*Run.NAtoms];
	hv = spd.x*spd.x + spd.y*spd.y + spd.z*spd.z;
	hv = hv * 1.0364269e2 * 0.5f * spd.w; //kinetic energy in eV!!!
	T[j*2+1] += hv;  
      }
      ClusterTemps[j*2+1] = 2.0f*T[j*2+1]/(3.0f*(ThermoBounds_h.w-ThermoBounds_h.z)) * 11604.505008098f;
    }
  }
  */

}


 
