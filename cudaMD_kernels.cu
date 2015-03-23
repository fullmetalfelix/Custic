// -*- c++ -*-





//MD functions
extern "C"
{

  // Kb = 1.3806503e-23 SI
  // 1 uma = 1.660538782e-27 Kg
  // 1 eV  = 1.60217646e-19 J
  // 1 eV/(Åuma) = 0.964853382e18  m/s² = 0.964853382e-2 Å/fs² (acceleration)
  // 1 N = 0.602214179e7 uma Å/fs²
  // 1 Å/fs = 

  //positions are in Å
  //speeds are in Å/fs
  //forces are in eV/Å
  //time is in fs
  //masses are in uma


  /*--- LEAPFROG INTEGRATION:
    update the velocities
    apply the thermostats to the velocities
    update the positions using the thermalized velocities
    copy the new positions to the gpu
   */
  

  __global__ void PosUpdate(float4 *X, float4 *V, float4 *H)
  {    
    float4 pos, spd;
       
    //get the right block index (separation between force and energy)
    int blockx = blockIdx.x;
    
    int imgidx = (blockx/BperIMG_d); //index of my image in shared to save space
    blockx -= imgidx*BperIMG_d;

    int idx = (blockx * BLOCKSIZE + threadIdx.x); //index of the atom this thread will take care of
    short isok = (idx < NAtoms_d);
    idx *= isok;
    int windex = idx + imgidx*NAtoms_d;
    
    pos = X[windex]; //position of the atom for which we will calculate the potential contribution
    spd = V[windex];
    
    //isok *= (H[windex].w >= 0);
    if(isok == 1){
      pos.x += MDstep_d * spd.x ;
      pos.y += MDstep_d * spd.y ;
      pos.z += MDstep_d * spd.z ;
      X[windex] = pos;
    }
  } 


  __global__ void SpdUpdate_plain(float3 *F, float4 *V, float4 *H)
  {    
    float4 spd;
    float3 frc;
    float imass;

    //get the right block index (separation between force and energy)
    int blockx = blockIdx.x;

    //get the image index
    int imgidx = (blockx/BperIMG_d); //index of my image in shared to save space
    blockx -= imgidx*BperIMG_d;

    int idx = (blockx * BLOCKSIZE + threadIdx.x); //index of the atom this thread will take care of
    short isok = (idx < NAtoms_d);
    idx *= isok;
    int windex = idx + imgidx*NAtoms_d;

    spd = V[windex];
    frc = F[windex];
    imass = 1.0f/spd.w;

    spd.x += MDstep_d * frc.x * 0.00964853382f * imass;
    spd.y += MDstep_d * frc.y * 0.00964853382f * imass;
    spd.z += MDstep_d * frc.z * 0.00964853382f * imass;

    if(H[windex].w < 0){
      spd.x = 0;
      spd.y = 0;
      spd.z = 0;
    }

    if(isok == 1){
      V[windex] = spd;
    }

  }


  __global__ void SpdUpdate(float3 *F, float4 *V, float4 *H, float *T)
  {    
    float4 spd;
    float3 frc;
    float lambda;
    float imass;

    //get the image index
    short imgidx = (blockIdx.x/BperIMG_d); //index of my image in shared to save space
    short idx = ((blockIdx.x-imgidx*BperIMG_d)*BLOCKSIZE + threadIdx.x); //index of the atom
    short isok = (idx < NAtoms_d);
    idx *= isok;
    short windex = idx + imgidx*NAtoms_d;
    
    //check if it is in the tip thermostat
    short thr = 0;
    if( (idx < ThermoBounds_d.y) && (idx >= ThermoBounds_d.x) )
      thr = 1;
    if( (idx < ThermoBounds_d.w) && (idx >= ThermoBounds_d.z) )
      thr = 2;
    
    thr--; //thr = 0 for tipthermo, 1 for surfthermo, -1 the rest
    
    //load the right temperature
    lambda = 1.0f;
    if(thr != -1){
      lambda += MDstep_d*Tcup_d[thr]*(Tset_d[thr]/T[imgidx*2+thr]-1.0f);
      //lambda = Tset_d[thr]/T[imgidx*2+thr];
    }
    lambda = sqrtf(lambda);
    
    
    spd = V[windex];
    frc = F[windex];
    imass = 0.00964853382f * MDstep_d /spd.w;
    
    //MDstep_d is in femtoseconds
    spd.x = spd.x*lambda + frc.x*imass;
    spd.y = spd.y*lambda + frc.y*imass;
    spd.z = spd.z*lambda + frc.z*imass;

    //if(thr != -1){
    //  spd.x = lambda;
    //  spd.y = Tset_d[thr];
    //  spd.z = T[imgidx*2+thr];
    //}

    if(H[windex].w < 0){
      spd.x = 0;
      spd.y = 0;
      spd.z = 0;
    }
  
    if(isok == 1){
      V[windex] = spd;
    }

  }

  /*
//DEBUG STUFF

  void SpdUpdate_cpu()
  {

    float lambda;
    float4 spd;
    float3 frc;
    float4 nspd[Run.NAtoms*Run.MDimages];

    cudaMemcpy(nspd,Speeds_d,f4_NObj,cudaMemcpyDeviceToHost);

    for(int j=0;j<Run.MDimages;j++)
      for(int i=0;i<Run.NAtoms;i++){
	lambda = 1.0f;
	frc = forces_h[j*Run.NAtoms+i];
	spd = Speeds_h[j*Run.NAtoms+i];

	if( (i<ThermoBounds_h.y) && (i>=ThermoBounds_h.x) )
	  lambda += Run.MD_Step/Run.TCoupTip*(Run.TempTip/Temperature_CPU[j*2]-1.0f);
	if( (i<ThermoBounds_h.w) && (i>=ThermoBounds_h.z) )
	  lambda += Run.MD_Step/Run.TCoupSurf*(Run.TempSurf/Temperature_CPU[j*2+1]-1.0f);
	  
	lambda *= 0.00964853382f * Run.MD_Step /spd.w;
	spd.x += frc.x * lambda;
	spd.y += frc.y * lambda;
	spd.z += frc.z * lambda;
	if(Holders_h[i+j*Run.NAtoms].w < 0){
	  spd.x = 0;
	  spd.y = 0;
	  spd.z = 0;
	}
	Speeds_h[j*Run.NAtoms+i] = spd;
      }

  }
  

  void SpeedUpdate()
  {

    unsigned int t1, t2;
    
    //backup copy of speeds
    float4 speedbk[Run.NAtoms*Run.MDimages];
    cudaMemcpy(speedbk,Speeds_d,f4_NObj,cudaMemcpyDeviceToHost);
    cudaMemcpy(Speeds_h,Speeds_d,f4_NObj,cudaMemcpyDeviceToHost);
    cudaMemcpy(Temperature_CPU,Temperature_GPU,2*Run.MDimages*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(forces_h,forces_d,f3_NObj,cudaMemcpyDeviceToHost);

    cutCreateTimer(&t1);
    cutStartTimer(t1);
    
    for(int i=0;i<1;i++)
      SpdUpdate_cpu( );

    cutStopTimer(t1);

    //save the result
    float4 speedscpu[Run.NAtoms*Run.MDimages];
    memcpy(speedscpu,Speeds_h,f4_NObj);
    
    
    //reset the gpu array
    cudaMemcpy(Speeds_d,speedbk,f4_NObj,cudaMemcpyHostToDevice);

    cutCreateTimer(&t2);
    cutStartTimer(t2);
    for(int i=0;i<1;i++){
      SpdUpdate<<<SlowGrid, SlowBlock>>>(forces_d, Speeds_d, Holders_d, Temperature_GPU);  //update the velocities
      cudaThreadSynchronize();
    }
    cutStopTimer(t2);

    printf("Speeds updated: %f - %f\n",cutGetTimerValue(t1)/100.0f,cutGetTimerValue(t2)/100.0f);

    //get the result
    cudaMemcpy(Speeds_h,Speeds_d,f4_NObj,cudaMemcpyDeviceToHost);
    
    FILE *fp = fopen("speeds.txt","w");
    for(int i=0;i<Run.NAtoms*Run.MDimages;i++){
      
      fprintf(fp,"%f %f %f\n",Speeds_h[i].x,Speeds_h[i].y,(Speeds_h[i].z-speedscpu[i].z));
    }
    fclose(fp);

    cutDeleteTimer(t1);
    cutDeleteTimer(t2);

  }

  */

}


 
