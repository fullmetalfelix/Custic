// -*- c++ -*-

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
//#include "defos.h"

extern "C"
{
  FILE* xyzMovie;
  FILE *MDforce;
  
  int elapsed;
  float elapsedt,time0;

  
  //compressed binary xyz movie functions
  void xyzMovie_Init(char* basename, int loop)
  {
    char filename[100];

    sprintf(filename,"%s.loop%i.ctj",basename,loop);
    xyzMovie = fopen(filename,"wb");
    fwrite(&Run.NAtoms,sizeof(int),1,xyzMovie);
    fwrite(&Run.MDimages,sizeof(int),1,xyzMovie);
  }

  /* *** Print a trajectory frame ***
  // trajectory file format is:
  // |Natoms|Nimages|
  // |int32 | int32 |
  //
  // |type| |x1 |y1 |z1 |Hx1|Hy1|Hz1| |x2 |y2 |z2 |Hx2|Hy2|Hz2| ... for each image
  // |5chr| |s16|s16|s16|s16|s16|s16| |s16|s16|s16|s16|s16|s16| ...
  // ... for each frame
  // ---------------------------- */
  void xyzMovie_printf( )
  {
    int i,j;  
    short3 pos;
    
    //fetch the positions and speeds from gpu memory
    cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost); //positions
    cudaMemcpy(Speeds_h,  Speeds_d,  f4_NObj, cudaMemcpyDeviceToHost); //speeds
    cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost); //holders positions
    
    
    //write the stuff in binary format
    for(i=0;i<Run.NAtoms;i++){
      fwrite(Types[System[i].TypeID].Name, sizeof(char),5, xyzMovie); //->typename as 5 chars
      
      for(j=0;j<Run.MDimages;j++){
	//compress the position
	pos.x = (short)ceil(Charges_h[i+j*Run.NAtoms].x*100.0f);
	pos.y = (short)ceil(Charges_h[i+j*Run.NAtoms].y*100.0f);
	pos.z = (short)ceil(Charges_h[i+j*Run.NAtoms].z*100.0f);
	//CAREFUL! HOLDERS HERE
	//pos.x = (short)ceil(Holders_h[i+j*Run.NAtoms].x*100.0f);
	//pos.y = (short)ceil(Holders_h[i+j*Run.NAtoms].y*100.0f);
	//pos.z = (short)ceil(Holders_h[i+j*Run.NAtoms].z*100.0f);

	fwrite(&pos,sizeof(short3),1,xyzMovie); //->position as short3
	//compress the holder position
	//pos.x = (short)ceil(Holders_h[i+j*Run.NAtoms].x*100.0f);
	//pos.y = (short)ceil(Holders_h[i+j*Run.NAtoms].y*100.0f);
	//pos.z = (short)ceil(Holders_h[i+j*Run.NAtoms].z*100.0f);
	//fwrite(&pos,sizeof(short3),1,xyzMovie); //->holders as short3!!!

	//write speed
	//fwrite(&Speeds_h[i+j*Run.NAtoms],sizeof(float3),1,xyzMovie); //->speed as float3!!!
	fwrite(&Holders_h[i+j*Run.NAtoms],sizeof(float3),1,xyzMovie); //->holders instead! as float3!!!

      } 
    }
    
  }


  void xyzMovie_Close()
  {
    fclose(xyzMovie);
  }

  //********************************************************************
  
  //MD Stats files
  void stats_Init(char* forbase, int loop)
  {    
    char filename[100];
      
    sprintf(filename,"%s.loop%i.out",forbase,loop);
    MDforce = fopen(filename,"wb");      

    elapsed = 0;
    elapsedt = 0;
    time0 = 0;
  }
  
  void stats_Close()
  {
    fclose(MDforce);  
  }




  //gathers values (E, Ek, T, ...) in the average collectors
  void StatsUpdate()
  {
    elapsed++;
    elapsedt += Run.MD_Step;
    
    GetTipForce(); //gpu calculation of tip forces

    //AvgTipPos.x += TipPos.x;
    //AvgTipPos.y += TipPos.y;
    //AvgTipPos.z += TipPos.z;

  }
  
  void ResetAverages( )
  {
    for(int j=0;j<Run.MDimages;j++){ //for each image
      TipForce_CPU[j] = make_float3(0.0f,0.0f,0.0f);
    }
    cudaMemcpy(TipForce_GPU, TipForce_CPU, Run.MDimages*sizeof(float3), cudaMemcpyHostToDevice );
    
    
    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos
    elapsed = 0;
    elapsedt= 0.0f;
  }
  

  /* *** Print the statistics ***
  // stats file was removed due to its uselessness
  // force file format is:
  // |time|xtip|ytip|ztip| |Fx1|Fy1|Fz1|T1t|T1s| |Fx2|Fy2|Fz2|T2t|T2s| ... |FxN|FyN|FzN|TNt|TNs|
  // |f32 |f32 |f32 |f32 | |f32|f32|f32|f32|f32| |f32|f32|f32|f32|f32| ... |f32|f32|f32|f32|f32|
  // ---------------------------- */
  void stats_printf(int step, int DebugPrint)
  {
    float time = time0+(step*Run.MD_Step - time0)*0.5f; //elapsedt/elapsed; //in ps
    
    if(elapsedt == 0.0f) //no printing if no time elapsed
      return;
    
    time0 = step*Run.MD_Step; //the start time is updated
    time *= 0.001f;  //time in ps
    
    //copy the temperatures and forces
    cudaMemcpy(Temperature_CPU, Temperature_GPU, 2*Run.MDimages*sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy(TipForce_CPU, TipForce_GPU, Run.MDimages*sizeof(float3), cudaMemcpyDeviceToHost );
    //----------------------------------

    //float4 zholder;
    cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost); //holders positions
    FirstHolder();
    cudaMemcpy(&(AvgTipPos.z), TipPosZ_d, sizeof(float), cudaMemcpyDeviceToHost);
    AvgTipPos.z /= elapsed;
    
    //print to file
    if(DebugPrint == true)
      printf("MD[0]: t = %8.5f, Tt=%8.5f, Ts=%8.5f, ztip=%6.2f fx=%6.3f fy=%6.3f fz=%6.3f | fy=%6.3f fz=%6.3f | fy=%6.3f fz=%6.3f |", time, 
	     Temperature_CPU[0], Temperature_CPU[1], AvgTipPos.z,
	     TipForce_CPU[0].x/elapsed, TipForce_CPU[0].y/elapsed, TipForce_CPU[0].z/elapsed
	     , TipForce_CPU[1].y/elapsed, TipForce_CPU[1].z/elapsed
	     , TipForce_CPU[2].y/elapsed, TipForce_CPU[2].z/elapsed );
    
    if(DebugPrint == true)
      printf("\n");
    

    //---------- force part --------------

    //print the shit
    fwrite(&time,sizeof(float),1,MDforce);         //-> time(averaged) as 1 float
    fwrite(&AvgTipPos,sizeof(float3),1,MDforce);   //-> tip position(averaged) as 1 float3

    //forces on the tip clusters only!
    for(int j=0;j<Run.MDimages;j++){
      
      //average the tip forces:
      TipForce_CPU[j].x /= elapsed;
      TipForce_CPU[j].y /= elapsed;
      TipForce_CPU[j].z /= elapsed;
      
      fwrite(&TipForce_CPU[j],sizeof(float3),1,MDforce); //->tip force(avg) as 1 float3
      fwrite(&Temperature_CPU[2*j],sizeof(float),1,MDforce); //->TIP temperature (not avg) as 1 float
      fwrite(&Temperature_CPU[2*j+1],sizeof(float),1,MDforce); //->SURF temperature (not avg) as 1 float

      //reset on cpu
      TipForce_CPU[j] = make_float3(0.0f,0.0f,0.0f);
      //Temperature_CPU[j] = 0.0f;
     }
    
    //reset the average positions
    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos
    
    //cudaMemcpy(Temperature_GPU, Temperature_CPU, Run.MDimages*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy(TipForce_GPU, TipForce_CPU, Run.MDimages*sizeof(float3), cudaMemcpyHostToDevice );
    
    elapsed = 0;
    elapsedt= 0.0f;
    cudaMemcpy(TipPosZ_d, &(elapsedt), sizeof(float), cudaMemcpyHostToDevice);

  }
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------


  //restart file management
  void WriteSimpleRestart(char* basename, int loop)
  {
    FILE *fp;
    char filename[100];
    
    sprintf(filename, "%s.loop%i.ctj",basename,loop);

    fp = fopen(filename, "wb");

    //copy back from gpu
    cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost );
    cudaMemcpy(Speeds_h,  Speeds_d,  f4_NObj, cudaMemcpyDeviceToHost );
    cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost );

    fwrite(Charges_h,sizeof(float4),Run.NAtoms*Run.MDimages,fp);
    fwrite(Holders_h,sizeof(float4),Run.NAtoms*Run.MDimages,fp);
    fwrite(Speeds_h, sizeof(float4),Run.NAtoms*Run.MDimages,fp);

    fclose(fp);

    printf("INFO! Restart file written.\n\n");
  }
  void ReadSimpleRestart(char* basename)
  {
    FILE *fp;
     
    //reload the given configuration
    fp = fopen(basename,"rb");

    fread(Charges_h,sizeof(float4),Run.NAtoms*Run.MDimages,fp);
    fread(Holders_h,sizeof(float4),Run.NAtoms*Run.MDimages,fp);
    fread(Speeds_h, sizeof(float4),Run.NAtoms*Run.MDimages,fp);
    fclose(fp);
    
    //send them to gpu
    cudaMemcpy(Charges_d, Charges_h, f4_NObj, cudaMemcpyHostToDevice);
    cudaMemcpy(Holders_d, Holders_h, f4_NObj, cudaMemcpyHostToDevice);
    cudaMemcpy(Speeds_d,  Speeds_h,  f4_NObj, cudaMemcpyHostToDevice);

    printf("INFO! Restart file (%s) loaded.\n\n",basename);
  }
  
  //*************************************************************************


  //VMD xyz format single frame outputs: one frame for each image
  void WriteCoords(char *basename)
  {
    FILE *fp;
    int i, img;
    char filename[100];
    
    sprintf(filename,"%s.xyz",basename);
    fp = fopen(filename,"w");
    
    //fetch from gpu
    cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost );
    cudaMemcpy(Speeds_h,  Speeds_d,  f4_NObj, cudaMemcpyDeviceToHost );
    //cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost );

    //print a vmd xyz file where each frame is the configuration of an image!

    fprintf(fp, "%i\n\n", Run.NAtoms);
    
    for(img=0;img<Run.MDimages;img++){
      for(i=0;i<Run.NAtoms;i++){
	
	fprintf(fp, "%s %f %f %f  %f %f %f\n", Types[System[i].TypeID].Name, 
		Charges_h[i+img*Run.NAtoms].x, Charges_h[i+img*Run.NAtoms].y, Charges_h[i+img*Run.NAtoms].z,
		Speeds_h[i+img*Run.NAtoms].x, Speeds_h[i+img*Run.NAtoms].y, Speeds_h[i+img*Run.NAtoms].z);
	
      }
      fprintf(fp,"\n\n");
    }
    fclose(fp);
    
  }


}
