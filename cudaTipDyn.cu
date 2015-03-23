// -*- c++ -*-


extern "C"
{

  void MoveCluster(float3 step, int ID);
  void ResetHolder( );


  void LinearMove( )
  {
    float3 dr;

    dr.x = Run.TipDynPars[3]*Run.TipDynPars[0]*Run.MD_Step*0.001f;
    dr.y = Run.TipDynPars[3]*Run.TipDynPars[1]*Run.MD_Step*0.001f;
    dr.z = Run.TipDynPars[3]*Run.TipDynPars[2]*Run.MD_Step*0.001f;
    
    //move the tip clusters - ALL OF THEM!
    MoveTipSmart( dr ); //move only the holders
    
  }
  
  float HiHolderZ, HiHolderZREF, HiHolderZDELTA;


  //this finds the coordinate of the first holder in the tip "holder" cluster
  //and assign it to TipHld
  //the holders need to be in the CPU memory
  void FirstHolder( )
  {
    
    TipHld.x = Holders_h[0].x;
    TipHld.y = Holders_h[0].y;
    TipHld.z = Holders_h[0].z;

  }



  //initialize the harmonic motion
  void HarmonicInit( )
  {
    float f = 2.0f*PI*Run.TipDynPars[0]; //in MHz - frequency is converted to angular pulse
    float period = 2.0f*PI/f;
    float acs, amp, thr;
    
    printf("Initializing harmonic motion...\n");

    thr = Run.TipDynPars[2];
    amp = Run.TipDynPars[1];

    //find the final time
    acs = acos(thr/amp - 1.0f)/f; //start time! in microseconds
    
    period = period -2*acs;
    period *= 1.0e9/Run.MD_Step; //in steps
    
    //printf("period %ffs = %isteps\n",period,(int)ceil(period));
    Run.MD_NSteps = (int)ceil(period); //put the right amount of steps to complete one cutted loop
    printf("  number of MD steps: %i - %f ns\n", Run.MD_NSteps, ceil(period*Run.MD_Step)*1.0e-6 );
    
    //store the corrected values in the tipdyn
    Run.TipDynPars[0] = f;
    Run.TipDynPars[2] = acs;
    
    
    
    //copy the tip dynamics parameters to gpu constant memory
    float3 tipdyn = make_float3( Run.TipDynPars[0],Run.TipDynPars[1],Run.TipDynPars[2] );
    cudaMemcpyToSymbol( TipHrm_d, &tipdyn, sizeof(float3), 0, cudaMemcpyHostToDevice);

    f = TipHldDist.z+Run.TipXY[2];
    cudaMemcpyToSymbol( TipOffset_d, &f, sizeof(float), 0, cudaMemcpyHostToDevice);

    f=0.0f;
    cudaMemcpy(TipPosZ_d, &f, sizeof(float), cudaMemcpyHostToDevice);

  }



  void HarmonicReInit( float3 actual, float3 ref )
  {
    float3 dr;

    //we should see where the tip holder effectively is
    
    //the holders are in CPU since at the end of the loop a frame is output
    //and and all data is copied from GPU memory!
    //at this point TipHld should contain the holder position after the loop
    
    //printf("The loop ended @ %f. REF %f\n",TipPos.z,Run.TipXY[2]);

    dr.x = 0.0f;
    dr.y = 0.0f;
    dr.z = ref.z - actual.z;
    
    //move the tip clusters (atoms and holders) on the GPU
    MoveTip( dr );

    //cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost);
    //cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost);

    //TipPos.z = 0.0f;//Run.TipXY[2];
    ResetHolder( ); //make sure the positions are uberprecise and copy to cpu
    cudaMemcpy(TipPosZ_d, &(dr.x), sizeof(float), cudaMemcpyHostToDevice);
  }

  //initialize the harmonic motion
  void HarmonicTRInit( )
  {
    float f = 2.0f*PI*Run.TipDynPars[0]; //in MHz - frequency is converted to angular pulse
    float period = 2.0f*PI/f;
    float acs, amp, thr;
    
    printf("Initializing harmonic motion...\n");

    thr = Run.TipDynPars[2];
    amp = Run.TipDynPars[1];

    //find the final time
    acs = acos(thr/amp - 1.0f)/f; //start time! in microseconds
    
    period = period -2*acs;
    period *= 1.0e9/Run.MD_Step; //in steps
    
    //printf("period %ffs = %isteps\n",period,(int)ceil(period));
    Run.MD_NSteps = (int)ceil(period); //put the right amount of steps to complete one cutted loop
    printf("  number of MD steps: %i - %f ns\n", Run.MD_NSteps, ceil(period*Run.MD_Step)*1.0e-6 );
    
    //store the corrected values in the tipdyn
    Run.TipDynPars[0] = f;   //the angular pulse
    Run.TipDynPars[2] = acs; //?
                             //Pars[1] is the amplitude in Å
    
    
    //copy the tip dynamics parameters to gpu constant memory
    
    float3 tipdyn = make_float3( Run.TipDynPars[0],Run.TipDynPars[1],Run.TipDynPars[2] );
    cudaMemcpyToSymbol( TipHrm_d, &tipdyn, sizeof(float3), 0, cudaMemcpyHostToDevice);

    f = TipHldDist.y+Run.TipXY[1];
    cudaMemcpyToSymbol( TipOffset_d, &f, sizeof(float), 0, cudaMemcpyHostToDevice);

    f=0.0f;
    cudaMemcpy(TipPosZ_d, &f, sizeof(float), cudaMemcpyHostToDevice);

  }
  void HarmonicTRReInit( float3 actual, float3 ref )
  {
    float3 dr;

    //we should see where the tip holder effectively is
    
    //the holders are in CPU since at the end of the loop a frame is output
    //and and all data is copied from GPU memory!
    //at this point TipHld should contain the holder position after the loop
    
    //printf("The loop ended @ %f. REF %f\n",TipPos.z,Run.TipXY[2]);

    dr.x = 0.0f;
    dr.z = 0.0f;
    dr.y = ref.y - actual.y;
    
    //move the tip clusters (atoms and holders) on the GPU
    MoveTip( dr );

    //cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost);
    //cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost);

    //TipPos.z = 0.0f;//Run.TipXY[2];
    ResetHolder( ); //make sure the positions are uberprecise and copy to cpu
    cudaMemcpy(TipPosZ_d, &(dr.x), sizeof(float), cudaMemcpyHostToDevice);
  }




  void TipDynamics(int step)
  {
    
    float time = Run.MD_Step*step*1.0e-9;
    
    if(Run.TipDyn == 0) //static
      return;

    if(Run.TipDyn == 1){ //linear
      //printf("ASD\n");
      LinearMove( );
      return;
    }
    
    if(Run.TipDyn == 2){ //harmonic
      //printf("ASD\n");
      //HarmonicMove(step);
      MoveTipHarmonic( time );
      return;
    }

    if(Run.TipDyn == 3){ //TR  ACHTUNG!
      MoveTipHarmonicTR_v2( time );
    }
    
  }



}
