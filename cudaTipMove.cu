// -*- c++ -*-


extern "C"
{
  //********************************************************************************
  /* Tip movement:
     this kernel moves the tip atoms and their holders.
   */
  __global__ void MoveTip1_GPU( float4 *X, float4 *H, short *isTip, float3 step)
  {
    float4 x,h;
 
    short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=istip
    asd.x = (blockIdx.x/FperIMG_d);
    asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.z = asd.y * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms
 
    short isok = (asd.z < NAtoms_d);
    asd.z *= isok;
  
    x = X[asd.z+asd.x];
    h = H[asd.z+asd.x];
    asd.w = isTip[asd.z];

    
    //move the holder
    if(asd.w == 1){
      h.x += step.x;//*asd.w;
      h.y += step.y;//*asd.w;
      h.z += step.z;//*asd.w;
      
      x.x += step.x;//*asd.w;
      x.y += step.y;//*asd.w;
      x.z += step.z;//*asd.w;
    }

    if(isok == 1){
      X[asd.z+asd.x] = x;
      H[asd.z+asd.x] = h;
    }  
  }
  
  //********************************************************************************
  //********************************************************************************
  /* Tip movement:
     this kernel moves the tip atoms and their holders.
   */
  __global__ void MoveTip_GPU( float4 *X, float4 *H, float3 step)
  {
    float4 x,h;
 
    short2 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=istip
    asd.x = (blockIdx.x/FperIMG_d);
    //asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.y = (blockIdx.x - asd.x*FperIMG_d) * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms
 
    short isok = (asd.y < TipBound_d);
    asd.y *= isok;
  
    x = X[asd.y+asd.x];
    h = H[asd.y+asd.x];

    
    //move the holder
    h.x += step.x;//*asd.w;
    h.y += step.y;//*asd.w;
    h.z += step.z;//*asd.w;
    
    x.x += step.x;//*asd.w;
    x.y += step.y;//*asd.w;
    x.z += step.z;//*asd.w;
  

    if(isok == 1){
      X[asd.y+asd.x] = x;
      H[asd.y+asd.x] = h;
    }  
  }
  void MoveTip( float3 step )
  {
    MoveTip_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, step);
    cutilCheckMsg("Kernel execution failed");
    cudaThreadSynchronize();
  }




  //********************************************************************************
  /* Smart tip movement kernel
     this kernel moves the tip atoms in a smart way by step:
     The atomic holders are moved, but the real atoms are moved only if
     they are fixed. This way atoms pinned by a spring will not be moved
     directly by this kernel.    
   */
  __global__ void MoveTipSmart1_GPU( float4 *X, float4 *H, short *isTip, float3 step)
  {
    float4 x,h;
    short tip;

    short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/FperIMG_d);
    asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.z = asd.y * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms

    short isok = (asd.z < NAtoms_d);
    asd.z *= isok;
  
    x = X[asd.z+asd.x];
    h = H[asd.z+asd.x];
    tip = isTip[asd.z];
    
    //the holder is moved anyway!
    h.x += step.x*tip;
    h.y += step.y*tip;
    h.z += step.z*tip;
    
    tip *= (h.w < 0.0f); //tip=1 if fixed, tip=0 if free

    //if the atom is fixed and is in the tip, it has to move with its holder
    x.x += step.x*tip;
    x.y += step.y*tip;
    x.z += step.z*tip;

    if(isok == 1){
      X[asd.z+asd.x] = x;
      H[asd.z+asd.x] = h;
    }  
  }
  __global__ void MoveTipSmart_GPU( float4 *X, float4 *H, float3 step)
  {
    float4 x,h;
    short tip;

    short2 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/FperIMG_d);
    //asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.y = (blockIdx.x - asd.x*FperIMG_d) * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms

    short isok = (asd.y < TipBound_d);
    asd.y *= isok;
  

    
    //the holder is moved anyway
    if(isok == 1){

      x = X[asd.y+asd.x];
      h = H[asd.y+asd.x];

      h.x += step.x;
      h.y += step.y;
      h.z += step.z;

      tip = (h.w < 0.0f); //tip=1 if fixed, tip=0 if free

      //if the atom is fixed and is in the tip, it has to move with its holder
      x.x += step.x*tip;
      x.y += step.y*tip;
      x.z += step.z*tip;

      X[asd.y+asd.x] = x;
      H[asd.y+asd.x] = h;
    }
    
  }

  //********************************************************************************
  /* Smart tip movement kernel 2 (Harmonic)
     this kernel moves the tip atoms in a smart way:
     The atomic holders are moved, but the real atoms are moved only if
     they are fixed. This way atoms pinned by a spring will not be moved
     directly by this kernel.
     The displacement is calculated by the first thread in each block from the value
     of the timestep "tstep" and the fixed distance "thdist" between the tip sticker
     and the first holder + ztip0.
   */
  __global__ void MoveTipHarmonic_GPU( float4 *X, float4 *H, float tt, float* ztip)
  {
    float4 x,h;
    short tip;
     __shared__ float step;
    
    short2 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/FperIMG_d);
    //asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.y = (blockIdx.x - asd.x*FperIMG_d) * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms
    
    short isok = (asd.y < TipBound_d);
    asd.y *= isok;
    
    //the first thread computes the displacement
    if(threadIdx.x == 0){
      step = -2.0f*TipHrm_d.y*(sin(TipHrm_d.x*(0.5f*tt+TipHrm_d.z))*sin(TipHrm_d.x*0.5f*tt));
      if(blockIdx.x == 0)
	ztip[0] += step;
      step += TipOffset_d; //this is where the tip holder should be
      step -= H[0].z; //this is how much the holder should be displaced
    }
    __syncthreads();
    
    if(isok == 1){
      
      x = X[asd.y+asd.x];
      h = H[asd.y+asd.x];
      
      //the holder is moved anyway!
      h.z += step;
      
      tip = (h.w < 0.0f); //tip=1 if fixed, tip=0 if free
    
      //if the atom is fixed and is in the tip, it has to move with its holder
      x.z += step*tip; //ONLY LATERAL MOTION
      
      X[asd.y+asd.x] = x;
      H[asd.y+asd.x] = h;
    }  
  }

  __global__ void MoveTipHarmonicTR_GPU( float4 *X, float4 *H, float tt, float* ytip)
  {
    float4 x,h;
    short tip;
     __shared__ float step;
    
    short2 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/FperIMG_d);
    //asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.y = (blockIdx.x - asd.x*FperIMG_d) * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms
    
    short isok = (asd.y < TipBound_d); //check if the thread
    asd.y *= isok;
    
    //the first thread computes the displacement
    if(threadIdx.x == 0){
      step = -2.0f*TipHrm_d.y*(sin(TipHrm_d.x*(0.5f*tt+TipHrm_d.z))*sin(TipHrm_d.x*0.5f*tt));
      if(blockIdx.x == 0)
	ytip[0] += step;
      step += TipOffset_d; //this is where the tip holder should be
      step -= H[0].y; //this is how much the holder should be displaced
    }
    __syncthreads();
    
    if(isok == 1){
      
      x = X[asd.y+asd.x];
      h = H[asd.y+asd.x];
      
      //the holder is moved anyway!
      h.y += step;
      
      tip = (h.w < 0.0f); //tip=1 if fixed, tip=0 if free
    
      //if the atom is fixed and is in the tip, it has to move with its holder
      x.y += step*tip; //ONLY VERTICAL MOTION
      
      X[asd.y+asd.x] = x;
      H[asd.y+asd.x] = h;
    }  
  }

  // --- lateral oscillation kernel -----------------------------------------------------------------
  __global__ void MoveTipHarmonicTR_GPU_v2( float4 *X, float4 *H, float tt, float* ytip, float4 *HR)
  { //mapping: 1 block per image, 1 thread per tip atom
    float4 x;
    float hw,hy;
    short tip;
    __shared__ float step,newH0y;
    int idx = threadIdx.x + blockIdx.x*NAtoms_d; //index of the atom to load
 
    //the first thread computes the displacement
    if( threadIdx.x == 0 ){

      x = H[blockIdx.x]; //position of the first holder
      hw= x.w;  //spring of the holder
      hy= x.y;  //y of the holder

      step = -2.0f*TipHrm_d.y*(sin(TipHrm_d.x*(0.5f*tt+TipHrm_d.z))*sin(TipHrm_d.x*0.5f*tt));
      ytip[blockIdx.x] += step;
      step += TipOffset_d; //this is where the tip holder should be
      step -= hy; //this is how much the (FIRST) holder should be displaced

      hy += step; //new position for the first holder
      
      newH0y = hy;
    }
    __syncthreads();
    
    hw=H[idx].w; //load the spring of the holder
    x=X[idx]; //load the atomic position
    
    hy = newH0y + HR[threadIdx.x].y;
    x.y += step*(hw < 0.0f); //move only if frozen!
    
    H[idx].y = hy;
    X[idx].y = x.y;    
  }

  void MoveTipHarmonicTR_v2( float time )
  {   
    MoveTipHarmonicTR_GPU_v2<<<Run.MDimages, Run.NTipAtoms>>>(Charges_d, Holders_d, time, TipPosZ_d, HolderR_d);
    cudaThreadSynchronize();
  }
  //------------------------------------------------------------------------------------------

  void MoveTipSmart( float3 step )
  {
    //MoveTipSmart_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, TipAtoms_d, step);
    MoveTipSmart_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, step);
    cudaThreadSynchronize();
  }

  void MoveTipHarmonic( float time )
  {
    MoveTipHarmonic_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, time, TipPosZ_d);
    cudaThreadSynchronize();
  }



  // --- lateral oscillation kernel -----------------------------------------------------------------
  __global__ void MoveTipHarmonicTR_BITS_GPU( float4 *X, float4 *H, //positions and holders
					      float tt, float* ytip, float4 *HR)
  { //mapping: 1 block per image, 1 thread per tip atom
    float4 x;
    float hw,hy;
    short tip;
    __shared__ float step,newH0y;
    int idx = threadIdx.x + blockIdx.x*NAtoms_d; //index of the atom to load
 
    //the first thread computes the displacement
    if( threadIdx.x == 0 ){

      x = H[blockIdx.x]; //position of the first holder
      hw= x.w;  //spring of the holder
      hy= x.y;  //y of the holder

      step = -2.0f*TipHrm_d.y*(sin(TipHrm_d.x*(0.5f*tt+TipHrm_d.z))*sin(TipHrm_d.x*0.5f*tt));
      ytip[blockIdx.x] += step;
      step += TipOffset_d; //this is where the tip holder should be
      step -= hy; //this is how much the (FIRST) holder should be displaced

      hy += step; //new position for the first holder
      
      newH0y = hy;
    }
    __syncthreads();
    
    hw=H[idx].w; //load the spring of the holder
    x=X[idx]; //load the atomic position
    
    hy = newH0y + HR[threadIdx.x].y;
    x.y += step*(hw < 0.0f); //move only if frozen!
    
    H[idx].y = hy;
    X[idx].y = x.y;    
  }

  //void MoveTipHarmonicTR_BITS_v2( float time )
  //  {   
  //    MoveTipHarmonicTR_BITS_GPU_v2<<<Run.MDimages, Run.NTipAtoms>>>(Charges_d, Holders_d, time, TipPosZ_d, HolderR_d);
  //    cudaThreadSynchronize();
  //  }
  //------------------------------------------------------------------------------------------













  //void MoveTipHarmonicTR( float time )
  //{
  //  MoveTipHarmonicTR_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, time, TipPosZ_d);
  //  cudaThreadSynchronize();
  //}
  











  /*
  

  void MoveTip_test( )
  {
    float time = 0.0f;
    float3 step = make_float3(1.0f,0.0f,0.0f);
    dim3 blok(BLOCKSIZE);
    dim3 grid(Run.BperIMG * Run.MDimages);
    if(grid.x == 0) grid.x = 1;

    unsigned int t1,t2;

    
    cudaMemcpy(ChargeT_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost); // back up in T


    cutCreateTimer(&t1);
    cutStartTimer(t1);

    //int i;
    for(int i=0;i<1000;i++){
      //MoveTip_GPU<<<grid, blok>>>(Charges_d, Holders_d, TipAtoms_d, step);
      //MoveTip_GPU_hld<<<grid, blok>>>( Holders_d, TipAtoms_d, step);
      MoveTipSmart_GPU<<<grid, blok>>>(Charges_d, Holders_d, step);
      cudaThreadSynchronize();
    }

    cutStopTimer(t1);
    //printf("moved: %f ms\n",cutGetTimerValue(t1)/1000);
    
    
    cudaMemcpy(ChargeO_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost);
    
    //do it on the cpu
    cudaMemcpy(Charges_d, ChargeT_h, f4_NObj, cudaMemcpyHostToDevice); // restore in gpu from T
    cutCreateTimer(&t2);
    cutStartTimer(t2);

    for(int i=0;i<1000;i++){
      MoveTipHarmonic_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, time, TipPosZ_d);
      cudaThreadSynchronize();
    }
    /*
    for(int loop=0;loop<1000;loop++)
      for(int j=0;j<Run.MDimages;j++)
	for(int i=0;i<Run.NAtoms;i++){
	  Charges_h[i+j*Run.NAtoms].x += step.x*TipAtoms_h[i];
	  Charges_h[i+j*Run.NAtoms].y += step.y*TipAtoms_h[i];
	  Charges_h[i+j*Run.NAtoms].z += step.z*TipAtoms_h[i];
	}
  ///* /
    cutStopTimer(t2);
    
    
    printf("moved: %f - %f ms\n",cutGetTimerValue(t2)/1000,cutGetTimerValue(t1)/1000);
    
    /*
    FILE *fp = fopen("poscheck.txt","w");

    for(int j=0;j<Run.MDimages;j++)
      for(int i=0;i<Run.NAtoms;i++){
	fprintf(fp,"%f %f %f\n",0*Charges_h[i+j*Run.NAtoms].x-ChargeO_h[i+j*Run.NAtoms].x,
		Charges_h[i+j*Run.NAtoms].y-ChargeO_h[i+j*Run.NAtoms].y,
		Charges_h[i+j*Run.NAtoms].z-ChargeO_h[i+j*Run.NAtoms].z);
      }

    fclose(fp);
    ///* /

    cutDeleteTimer(t1);
    cutDeleteTimer(t2);

  }
  


  //harmonic motion test
  



  /*
  void HarmonicTest( )
  {
    float time, tip, disp,z;
    FILE *fp = fopen("tiptraj","w");
    unsigned int t1,t2;

    //cutCreateTimer(&t1);
    //cutStartTimer(t1);

    z = 20.0f;
    fprintf(fp,"%f %f %f\n",0.0f,tip,z);
    for(int t=1;t<5000000;t++){
      time = Run.MD_Step*t*1.0e-9; //time in us
      tip  = Run.TipDynPars[1]*( cos(Run.TipDynPars[0]*(time+Run.TipDynPars[2]))+1.0f ); //where the tip should be
      tip -= Run.TipDynPars[1]*( cos(Run.TipDynPars[0]*(Run.TipDynPars[2]))+1.0f );

      tip = -2.0f*Run.TipDynPars[1]*(sin(Run.TipDynPars[0]*(0.5f*time+Run.TipDynPars[2]))*sin(Run.TipDynPars[0]*0.5f*time));
      tip += 5.0f;
      //disp = tip - Run.TipDynPars[1]*( cos(Run.TipDynPars[0]*(Run.TipDynPars[2]))+1.0f ); //the displacement
      disp = tip-z;
      z += disp;
      
      if(t%1000 == 0){
	fprintf(fp,"%f %f %f\n",(float)t/1000000.0f,tip,z);
      }
      
    }
    fclose(fp);

    dim3 blok(BLOCKSIZE);
    dim3 grid(Run.BperIMG * Run.MDimages);
    if(grid.x == 0) grid.x = 1;

    cutCreateTimer(&t2);
    cutStartTimer(t2);

    //fp = fopen("tiptrajgpu","w");
    z = 10.0f;
    for(int t=1;t<5000;t++){
      time = Run.MD_Step*t*1.0e-9; //time in us
      tip  = Run.TipDynPars[1]*( cos(Run.TipDynPars[0]*(time+Run.TipDynPars[2]))+1.0f ); //where the tip should be
      disp = tip-z;

      //Harmony_GPU<<<grid, blok>>>(Charges_d, Holders_d, TipAtoms_d, time, z);
      cudaThreadSynchronize();

      z += disp;

      //if(t%1000 == 0){
      //	fprintf(fp,"%f %f %f\n",(float)t/1000000.0f,tip,z);
      //}
      
    }
    //fclose(fp);

    cutStopTimer(t2);
    printf("Harmtest gpu: %f \n",cutGetTimerValue(t2)/5000.0f);
    cutDeleteTimer(t2);

  }
  */



  // OLD KERNELS
  /*
  __global__ void MoveTip_GPU( float4 *X, float4 *H, short *isTip, float3 step)
  {
    float4 x,h;
    short tip;

    short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/BperIMG_d);
    asd.y = blockIdx.x - asd.x*BperIMG_d;
    asd.z = asd.y * BLOCKSIZE + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms

    short isok = (asd.z < NAtoms_d);
    asd.z *= isok;
  
    
    x = X[asd.z+asd.x];
    h = H[asd.z+asd.x];
    tip = isTip[asd.z];

    x.x += step.x*tip;
    x.y += step.y*tip;
    x.z += step.z*tip;

    if(isok == 1){
      X[asd.z+asd.x] = x;

    }  
  }
  
  //move the holders of tip atoms... the actual atoms will not be displaced!
  __global__ void MoveTip_GPU_hld( float4 *H, short *isTip, float3 step)
  {
    float4 h;
    short tip;

    short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=mytype
    asd.x = (blockIdx.x/BperIMG_d);
    asd.y = blockIdx.x - asd.x*BperIMG_d;
    asd.z = asd.y * BLOCKSIZE + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms

    short isok = (asd.z < NAtoms_d);
    asd.z *= isok;
  
    h = H[asd.z+asd.x];
    tip = isTip[asd.z];

    h.x += step.x*tip;
    h.y += step.y*tip;
    h.z += step.z*tip;

    if(isok == 1){
      H[asd.z+asd.x] = h;

    }  
  }
  */

}
