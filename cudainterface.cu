// -*- c++ -*-

#include "cudef.ch"

#include <stdio.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cutil_inline.h>


#include "cudaAlloc.cu"
#include "cudaVectors.cu"

//#include <uberkernel_energy.cu>
#include "uberkernel_forces2.cu"
//#include <old_forces2.cu>
//#include <uberkernel_both.cu>

//#include <cudaForce.cu>

//#include <cudaEnergy.cu> //for kinetic energy calculation only!

#include "cudaTipMove.cu"
#include "cudaCoordTools.cu"
#include "cudaTipDyn.cu"

#include "cudaT.cu"
#include "cudaF.cu"

#include "cudaMD_kernels.cu"
#include "cudaMD.cu"


#include "output_binary.cu"



extern "C"
{

  int MDrun( int );
  int MDrun_harmonic( int );
  int MDequilibrate ( int );


  //double FirstHolder( );

  void benchmark( void );


  //check if a GPU device is ok
  int ProbeDevice(int device)
  {
  
    int c;
    cudaError_t err;
    cudaDeviceProp prop;
  
    cudaGetDeviceCount(&c);

    if(c == 0)
      {
	printf("FATAL! No CUDA enabled device was found!\n");
	return false;
      }

    err = cudaGetDeviceProperties(&prop,device);
    if (err != 0)
      {
	printf("FATAL! cudaSetDevice ERROR %i \n", err);
	return false;
      }
    else
      {
	printf("(%s)\n",prop.name);
      }

    MyProps = prop;
    cutilCheckMsg("probe device failed");
    return 1;

  }


  //Quasi static force calculation on a grid
  int GridQSF()
  {
    //save the initial positions AS THEY ARE IN THE INPUT FILE!!!
    memcpy(ChargeINIT, Charges_h, f4_NObj);  


    //TODO...

    return true;
  }


  //where it all starts...
  int MainProcess()
  {
    srand((unsigned int)time(NULL));
    
    //write the input coordinates for vmd
    WriteCoords("CRD_INPUT");
    
    //if needed, reload the starting configuration from file
    if(Run.MD_reload == true){
      ReadSimpleRestart(Run.MD_reloadfile); //read the restart file
    }
    else{
      RandomizeSpeed(-1); //randomize the speeds
      
    }
    
    TipReposition(true); //put the tip in the right spot!
    

    //clear the averages
    //ResetAverages( );

    //initialize harmonic motion
    if(Run.TipDyn == 2) 
      HarmonicInit( );
    if(Run.TipDyn == 3) 
      HarmonicTRInit( );   

    //save the initial positions AFTER REPOSITION/reload!!!
    memcpy(ChargeINIT, Charges_h, f4_NObj);  


    //write initial coordinates AFTER REPOSITION/RELOAD!!!
    WriteCoords("CRD_INITS");

    
    FILE *fp = fopen("changes.log","w");fclose(fp); //reset the chain log file


    //****************************************************
    
    GetSystemTemp();
    //SpeedUpdate();
    //MoveTip_test( );
    //return true;

    float3 HolderREF;
    //*** now run the MD ***
    if(Run.TipDyn >= 2){    //HARMONIC TIP
      
      FirstHolder( ); HolderREF = TipHld; //save the starting holder point for reference
      

      for(int loop=1;loop<=(int)ceil(Run.TipDynPars[3]);loop++){

	
	if(Run.MD_reinit == false){
	  //normal operation: consecutive loops
	  
	  if( Run.MD_equiloop == true )
	    MDequilibrate(  loop+Run.MD_reload_num );  //run equilibration every loop
	  else
	    if( loop == 1 )
	      MDequilibrate(  loop+Run.MD_reload_num );  //run one equilibration
	  
	  MDrun_harmonic( loop+Run.MD_reload_num );  //run the main loop
	  if(Run.TipDyn == 2)
	    HarmonicReInit( TipHld, HolderREF );
	  if(Run.TipDyn == 3)
	    HarmonicTRReInit( TipHld, HolderREF );
	  
	}else{
	  //reinit operation: reinitialize the tip every loop

	  cudaMemcpy(Charges_d,ChargeINIT, f4_NObj, cudaMemcpyHostToDevice); //restore the tip
	  RandomizeSpeed(-1); //randomize the speeds

	  MDequilibrate(  loop+Run.MD_reload_num );  //run equilibration every loop

	  MDrun_harmonic( loop+Run.MD_reload_num );  //run the main loop
	  if(Run.TipDyn == 2)
	    HarmonicReInit( TipHld, HolderREF );
	  if(Run.TipDyn == 3)
	    HarmonicTRReInit( TipHld, HolderREF );


	}

      }
    }
    else{          //OTHER TIP DYNAMICS
      MDequilibrate( 0 ); //run one equilibration
      //MDrun( 0 );  //run the main loop
    }
    
    
    WriteSimpleRestart("RST.final", 0 );
    WriteCoords("CRD_FINAL");
    return true;
    
  }
  
  

  int MDequilibrate ( int loop )
  {
    printf("Equilibrating...\n ");

    ResetAverages();
    stats_Init("force.EQ",0);
    if(Run.MD_SaveXYZ == true){
      xyzMovie_Init("EQout",loop);
      xyzMovie_printf( );
    }

    GetSystemTemp();
    GetForcesEnergy( ); //get the forces for MD
    for(int i=1; i<=Run.MD_ESteps; i++){  //main loop
      
      LeapFrog( );     //integrate trajectory
      GetForcesEnergy();    //get new forces
      StatsUpdate();        //update statistics
      
      if(i%Run.MD_Stats == 0){ //save stats & forces
	stats_printf(i, true); //do not print to console  
      }
      
      if((i%Run.MD_XYZfreq == 0) && (Run.MD_SaveXYZ == true)) //save movie
	xyzMovie_printf( );

    }
    
    stats_Close();
    if(Run.MD_SaveXYZ == true){
      xyzMovie_printf( );
      xyzMovie_Close();
    }
    
    printf("Equilibration finished.\n");
    WriteSimpleRestart("RST.EQ", loop ); //copy back and write restart

    return true;
  }

  
  int MDrun_harmonic ( int loop )
  {

    int i;
    float tipz0;
 
    printf("Running MD (loop %i)...\n",loop);

    //save the initial position of the holder
    FirstHolder( );
    tipz0 = TipHld.z; //make a backup
    printf("The first holder is: %f %f %f\n",TipHld.x,TipHld.y,TipHld.z);
    printf("The tip is         : %f %f %f\n",TipPos.x,TipPos.y,TipPos.z);
    printf("DeltaZ is          : %f\n",TipHld.z-TipPos.z);

    //make a backup copy of the system
    memcpy(ChargeO_h, Charges_h, f4_NObj);
    ResetAverages( ); //reset the averages

    //initialize the output
    stats_Init("force.MD",loop);
    if(Run.MD_SaveXYZ == true){
      xyzMovie_Init("MDout",loop);
      xyzMovie_printf( );
    }
    
    //main MD loop
    GetSystemTemp();
    GetForcesEnergy( );   //get the forces for MD
    i=1;
    do{
      LeapFrog( );        //update speeds and positions
      TipDynamics(i);
      GetForcesEnergy( ); //get the forces again
      StatsUpdate();
	
      if(i%Run.MD_Stats == 0){ //save stats & forces
	stats_printf(i,true);
      }
      
      if((i%Run.MD_XYZfreq == 0) && (Run.MD_SaveXYZ == true)) //save movie
	xyzMovie_printf( );
	
      i++;
    }while(i<=Run.MD_NSteps);
    //while(TipPos.z < Run.TipXY[2]);

    stats_printf(i,true);

    //close the output files
    if(Run.MD_SaveXYZ == true){
      xyzMovie_printf( );
      xyzMovie_Close();
    }
    stats_Close();
    
    //copy back and save restart
    WriteSimpleRestart("RST.MD", loop);

    FirstHolder( );
    printf("Loop ended: Zini %f, Zf %f \n",tipz0,TipHld.z);

    
    if(Run.ChainReset == 1)
      CheckChain( loop );

    return true;
    
  } 


  //-----------------------------------------------------

  int MDrun ( int loop )
  {

    int i;
      
    printf("Running MD (loop %i)...\n",loop);

    //ResetAverages( ); //reset the averages
    
    
    //initialize the output
    stats_Init("MDforce",loop);
    if(Run.MD_SaveXYZ == true)
      {
	xyzMovie_Init("MDout",loop);
	xyzMovie_printf( );
      }
    
    //main MD loop
    GetForcesEnergy( );            //get the forces for MD
    for(i=1;i<=Run.MD_NSteps;i++)
      {
	LeapFrog( );            //update speeds and positions
	TipDynamics(i);
	GetForcesEnergy( );     //get the forces again
	StatsUpdate();
	
	if(i%Run.MD_Stats == 0) //save stats & forces
	  {
	    stats_printf(i,true);
	    //WriteRestart(Charges_h, Speeds_h, i);
	  }
	    
	if((i%Run.MD_XYZfreq == 0) && (Run.MD_SaveXYZ == true)) //save movie
	  xyzMovie_printf( );

      }
    stats_printf(i,true);
    //close the output files
    if(Run.MD_SaveXYZ == true)
      {
	xyzMovie_printf( );
	xyzMovie_Close();
      }
    stats_Close();
    
    if(Run.ChainReset == 1)
      CheckChain( loop );

    //printf("first Holder atoms: %f\n",FirstHolder());

    return true;
    
  }
  



  void benchmark()
  {

    unsigned int t1,t2;

    
    cutCreateTimer(&t1);
    cutCreateTimer(&t2);

    
    cutStartTimer(t1);
    //GetForces( true, true );
    //GetEnergy( true );
    cutStopTimer(t1);


    cutStartTimer(t2);
    GetForcesEnergy( );
    cutStopTimer(t2);


    printf("F+E in %f ms\n",cutGetTimerValue(t1));
    printf("FandE in %f ms\n",cutGetTimerValue(t2));


    cutDeleteTimer(t1);
    cutDeleteTimer(t2);


  }
  





}
