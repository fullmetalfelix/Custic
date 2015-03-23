// -*- c++ -*-

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
//#include "defos.h"

extern "C"
{
  FILE** xyzMovie;
  FILE** MDstats, **MDforce;
  //FILE** xyzMovieG;
  int elapsed;

  //  void xyzMovie_Init(char* basename);
  //  void xyzMovie_printf(int Rejected, float4 *Charges, float3 *Speeds);
  //  void xyzMovie_Close();

  float GetKinetic(int);
  float GetEnergy(int DebugPrint);
  void  GetClustersForce();


  void OutputAlloc()
  {
    xyzMovie = (FILE**) malloc(Run.MDimages*sizeof(FILE*));
    MDstats  = (FILE**) malloc(Run.MDimages*sizeof(FILE*));
    MDforce  = (FILE**) malloc(Run.MDimages*sizeof(FILE*));
    printf("Output files allocated.\n");
  }
  void OutputDealloc()
  {
    free(xyzMovie);
    free(MDstats);
    free(MDforce);
    printf("Output files deallocated.\n");
  }


 //VMD xyz movie functions
  void xyzMovie_Init(char* basename, int loop)
  {
    int j;
    char filename[100];

    for(j=0;j<Run.MDimages;j++){
      sprintf(filename,"%s.i%i.loop%i.xyz",basename,j,loop);
      xyzMovie[j] = fopen(filename,"w");
      fprintf(xyzMovie[j], "%i\n\n", Run.NAtoms);
    }
     
  }

  //writes a frame in the file
  void xyzMovie_printf(float4 *Charges, float3 *Speeds)
  {
    int i,j;  
    
    for(j=0;j<Run.MDimages;j++){
      for(i=0;i<Run.NAtoms;i++)
	{
	  fprintf(xyzMovie[j], "%s %f %f %f  %f %f %f\n",Types[System[i].TypeID].Name, 
		  Charges[i+j*Run.NAtoms].x, Charges[i+j*Run.NAtoms].y, Charges[i+j*Run.NAtoms].z,
		  Speeds[i+j*Run.NAtoms].x, Speeds[i+j*Run.NAtoms].y, Speeds[i+j*Run.NAtoms].z);
	  
	}
      fprintf(xyzMovie[j],"\n\n");
    }
  

  }


  void xyzMovie_Close()
  {
    int j;
    for(j=0;j<Run.MDimages;j++){
      fclose(xyzMovie[j]);
    }
  }


  
  //MD Stats files
  void stats_Init(char* stsbase, char* forbase, int loop)
  {    
    int j;
    char filename[100];
    for(j=0;j<Run.MDimages;j++)
      {
	sprintf(filename,"%s.i%i.loop%i.out",stsbase,j,loop);
	MDstats[j] = fopen(filename,"w");

	sprintf(filename,"%s.i%i.loop%i.out",forbase,j,loop);
	MDforce[j] = fopen(filename,"w");      

      }
    elapsed = 0;

  }
  
  void stats_Close()
  {
    int j;
    for(j=0;j<Run.MDimages;j++)
      {
	fclose(MDstats[j]);
	fclose(MDforce[j]);  
      }
  }




  //gathers values (E, Ek, T, ...) in the average collectors
  void StatsUpdate()
  {
    int i,j;

    elapsed++;
 
    for(j=0;j<Run.MDimages;j++){
      AvgE[j] += Etot[j]; //GetEnergy(false); // update the potential energy
      AvgEK[j]+= GetKinetic(j);     // update the kinetic energy
    }

    GetClustersForce();       //this updates the average cluster force
    for(j=0;j<Run.MDimages;j++)
      for(i=0;i<Run.NCls;i++)   //updates the clusters average temperature
	ClusterTempsAvg[i+j*Run.NCls] += ClusterTemps[i+j*Run.NCls]; //clustertemps are updated in the leapfrog!

    AvgTipPos.x += TipPos.x;
    AvgTipPos.y += TipPos.y;
    AvgTipPos.z += TipPos.z;


  }
  void ResetAverages( void )
  {
    int i,j;

    for(j=0;j<Run.MDimages;j++){
      AvgE[j] = 0.0f;
      AvgEK[j]= 0.0f;
 
      for(i=0;i<Run.NCls;i++)
	{
	  ClusterTempsAvg[i+j*Run.NCls] = 0.0f;
	  ClusterForce[i+j*Run.NCls]    = make_float3( 0.0f, 0.0f, 0.0f );
	}
    }

    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos
    elapsed = 0;
  }


  //print the statistics (energies and temperatures)
  void stats_printf(int step, int DebugPrint)
  {

    float time = (step-Run.MD_Stats/2)/1000.0f; //in ps
    int i,j;

    for(j=0;j<Run.MDimages;j++){
      AvgE[j] /= elapsed;
      AvgEK[j]/= elapsed;
    }

    //print to file
    if(DebugPrint == true)
      printf("MD[0]: t = %f, E=%f, Ek=%f, Etot=%f, ",time, AvgE[0], AvgEK[0], AvgE[0]+AvgEK[0] );

    for(j=0;j<Run.MDimages;j++)
      {
	fprintf(MDstats[j], "%f %f %f %f ",time, AvgE[j], AvgEK[j], AvgE[j]+AvgEK[j]);
	for(i=0;i<Run.NCls;i++)
	  {
	    if(Clusters[i].Thermostat != 0)
	      {
		ClusterTempsAvg[i+j*Run.NCls] /= elapsed;
		if(DebugPrint == true && j==1)
		  printf("T[%s] %f, ",Clusters[i].Name,ClusterTempsAvg[i+j*Run.NCls]);
		fprintf(MDstats[j], "%f ",ClusterTempsAvg[i+j*Run.NCls]);
	      }
	    ClusterTempsAvg[i+j*Run.NCls] = 0.0f;  //reset the average clstemp
	  }
	fprintf(MDstats[j],"\n");

	//reset the avg E
	AvgE[j] = 0;
	AvgEK[j]= 0;
      }
    
    if(DebugPrint == true)
      printf("\n");
 

    //---------- force part --------------
    AvgTipPos.x /= elapsed;
    AvgTipPos.y /= elapsed;
    AvgTipPos.z /= elapsed;
    AvgTipHld.x /= elapsed;
    AvgTipHld.y /= elapsed;
    AvgTipHld.z /= elapsed;

    //print the shit
    for(j=0;j<Run.MDimages;j++)
      {
	fprintf(MDforce[j], "%f %f %f %f ",time, AvgTipPos.x, AvgTipPos.y, AvgTipPos.z);
	
	for(i=1;i<Run.NCls;i++) // skip the default
	  {
	    
	    ClusterForce[i+j*Run.NCls] = Scalar_f3_f( ClusterForce[i+j*Run.NCls], 1.0f/elapsed );  //divide the average
	    
	    if(Clusters[i].Thermostat == 0) //prints the forces on the clusters with no thermostat
	      fprintf(MDforce[j]," %f %f %f ",ClusterForce[i+j*Run.NCls].x, 
		      ClusterForce[i+j*Run.NCls].y, ClusterForce[i+j*Run.NCls].z);    
	    
	    ClusterForce[i+j*Run.NCls] = make_float3( 0.0f, 0.0f, 0.0f ); //reset the average
	  }
	fprintf(MDforce[j],"\n");
      }
    
    //reset the average positions
    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos  

    elapsed = 0;

  }

    
  
  
  void WriteCoords(char *basename, float4 *Charges, float3 *Speeds)
  {

    
    FILE *fp;
    int i, img;
    char filename[100];
    
    for(img=0;img<Run.MDimages;img++)
      {
	sprintf(filename,"%s.i%i.xyz",basename,img);
	//printf("filename = %s\n",filename);
	fp = fopen(filename,"w");

	/*
	if(xyzFormat == false)
	  {
	    for(j=0;j<Run.NCls;j++)
	      {
		//fprintf(fp,"%s\n",Clusters[j].Name);
		for(i=0;i<Run.NAtoms;i++)
		  {
		    
		    if(System[i].ClusterID == j)
		      {
			fprintf(fp,"%s   %f %f %f\n",Types[System[i].TypeID].Name, 
				Charges[i+img*Run.NAtoms].x, Charges[i+img*Run.NAtoms].y, Charges[i+img*Run.NAtoms].z);
		      }
		    
		  }
	      }
	  }
	  else //else write it in XYZ format for vmd
	  {*/
	    
	fprintf(fp, "%i\n\n", Run.NAtoms);
	for(i=0;i<Run.NAtoms;i++)
	  {
	    
	    fprintf(fp, "%s %f %f %f  %f %f %f\n", Types[System[i].TypeID].Name, 
		    Charges[i+img*Run.NAtoms].x, Charges[i+img*Run.NAtoms].y, Charges[i+img*Run.NAtoms].z,
		    Speeds[i+img*Run.NAtoms].x, Speeds[i+img*Run.NAtoms].y, Speeds[i+img*Run.NAtoms].z);
	    
	    
	  }
	fprintf(fp,"\n\n");

	fclose(fp);
      }
    

	
  }
  




  void WriteRestart(float4 *X, float3 *V, int step, int cycle)
  {
    
    FILE* fp;
    int i, j;
    char filename[100];

    for(j=0;j<Run.MDimages;j++) //for each image
      {
	sprintf(filename,"CHECKPOINT.i%i.RST",j);
	fp = fopen(filename,"w");


	fprintf(fp,"time  %i\n", step);
	fprintf(fp,"cycle %i\n", cycle);
	
	fprintf(fp,"tip %f %f %f\n", TipPos.x, TipPos.y, TipPos.z);
	fprintf(fp,"time  %i\n", step);



	//write all atom coords and speeds
	for(i=0;i<Run.NAtoms;i++){

	}

	fclose(fp);

      }
     

    

  }



}
