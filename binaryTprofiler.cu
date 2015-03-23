// -*- c++ -*-

#include <stdio.h>
//#include <stdlib.h>
//#include <strings.h>
#define true 1
#define false 0


//extract a given loop of a given image
int StartsWith(char strIn[], char strRef[]);



int main (int argc, char* argv[])
{

  char filename[100], nameout[100],aname[5];
  int i,j,loop,img, natoms, nimages;
  FILE *fin, *fout;
  
  if(argc<4){
    printf("Insufficient parameters!\n");
    return 0;
  }
  
  
  loop = atoi(argv[2]);
  img  = atoi(argv[3]);

  printf("Extracting loop %i for image %i\n",loop,img);

  sprintf(nameout , "%s.loop%i.img%i.xyz",argv[1],loop,img);
  sprintf(filename,"%s.loop%i.ctj",argv[1],loop);
  printf("Opening %s\n",filename);
  fin = fopen(filename,"rb");
  fout= fopen(nameout ,"w");
  
  //read the number of atoms
  fread(&natoms,sizeof(int),1,fin);
  fprintf(fout,"%i\n\n",natoms);
  fread(&nimages,sizeof(int),1,fin);

  short3 pos;
  float3 speed;


  while(!feof(fin)){
    for(i=0;i<natoms;i++){ //loop all the atoms
      
      //read the name
      if((int)fread(aname,sizeof(char),5,fin) == 0)
	break;
      
      fprintf(fout,"%s ",aname);
 
      for(j=0;j<nimages;j++){
      
	fread(&pos,sizeof(short3),1,fin);  //->position as float3
	
	if(j==img-1){
	  speed.x = ((float)pos.x)/100.0f;
	  speed.y = ((float)pos.y)/100.0f;
	  speed.z = ((float)pos.z)/100.0f;
	  fprintf(fout,"%f %f %f\n",speed.x,speed.y,speed.z);
	}
	
	fread(&speed,sizeof(float3),1,fin);  //->speed as float3
	//if(j==img-1){
	//  fprintf(fout,"%f %f %f\n",speed.x,speed.y,speed.z);
	//}
	
      }
    }//all atoms written - for one frame
    fprintf(fout,"\n\n");
    //printf("frame done!\n");
  }
  

  fclose(fin);fclose(fout);





}





int StartsWith(char strIn[], char strRef[])
{
  int i;
  
  for(i=0;i<strlen(strRef);i++)
    {
      if(strIn[i]!=strRef[i])
	return false;
      //printf("%c %c\n",strIn[i],strRef[i]);
    }
  return true;
}








/*

 //VMD xyz movie functions
  void xyzMovie_Init(char* basename, int loop)
  {
    char filename[100];

    sprintf(filename,"%s.loop%i.ctj",basename,loop);
    xyzMovie = fopen(filename,"wb");
    fwrite(&Run.NAtoms,sizeof(int),1,xyzMovie);
    
  }

  //writes a frame in the file
  void xyzMovie_printf(float4 *Charges, float3 *Speeds)
  {
    int i,j;  
    short3 pos;
    

    for(i=0;i<Run.NAtoms;i++){
      fwrite(Types[System[i].TypeID].Name, sizeof(char),5, xyzMovie); //->typename as 5 chars
      
      for(j=0;j<Run.MDimages;j++){
	//compress the position
	pos.x = (short)ceil(Charges[i+j*Run.NAtoms].x*100.0f);
	pos.y = (short)ceil(Charges[i+j*Run.NAtoms].y*100.0f);
	pos.z = (short)ceil(Charges[i+j*Run.NAtoms].z*100.0f);
	fwrite(&pos,sizeof(short3),1,xyzMovie);  //->position as float3

	fwrite(&Speeds[i+j*Run.NAtoms],sizeof(float3),1,xyzMovie); //->speed as float3
	
	//fwrite(&forces_h[i+j*Run.NAtoms],sizeof(float3),1,xyzMovie); //->force as float3

      } 
    }
    
  }


  void xyzMovie_Close()
  {
    fclose(xyzMovie);
  }


  
  //MD Stats files
  void stats_Init(char* stsbase, char* forbase, int loop)
  {    
    char filename[100];
    
    sprintf(filename,"%s.loop%i.out",stsbase,loop);
    MDstats = fopen(filename,"wb");
    
    sprintf(filename,"%s.loop%i.out",forbase,loop);
    MDforce = fopen(filename,"wb");      

    elapsed = 0;
    elapsedt = 0;
    time0 = 0;
  }
  
  void stats_Close()
  {
    fclose(MDstats);
    fclose(MDforce);  
  }




  //gathers values (E, Ek, T, ...) in the average collectors
  void StatsUpdate()
  {
    int i,j;

    elapsed++;
    elapsedt += Run.MD_Step;

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

    for(j=0;j<Run.MDimages;j++){ //for each image
      AvgE[j] = 0.0f;
      AvgEK[j]= 0.0f;
 
      for(i=0;i<Run.NCls;i++){ //loop over the clusters
	  ClusterTempsAvg[i+j*Run.NCls] = 0.0f;
	  ClusterForce[i+j*Run.NCls]    = make_float3( 0.0f, 0.0f, 0.0f );
	}
    }

    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos
    elapsed = 0;
    elapsedt= 0.0f;
  }


  //print the statistics (energies and temperatures)
  void stats_printf(int step, int DebugPrint)
  {
    float time = time0+(step*Run.MD_Step - time0)*0.5f; //elapsedt/elapsed; //in ps
    int i,j;

    time0 = step*Run.MD_Step; //the start time is updated
    time *= 0.001f;  //time in ps

    for(j=0;j<Run.MDimages;j++){
      AvgE[j] /= elapsed;
      AvgEK[j]/= elapsed;
    }

    //print to file
    if(DebugPrint == true)
      printf("MD[0]: t = %8.5f, E=%8.5f, Ek=%8.5f, ",time, AvgE[0], AvgEK[0], elapsedt );

    fwrite(&time,sizeof(float),1,MDstats);         //-> time(averaged) as 1 float
    
    for(j=0;j<Run.MDimages;j++){  //for each image
 
      for(i=0;i<Run.NCls;i++){  //for each cluster
	
	if(Clusters[i].Thermostat != 0){  //... with a thermostat/meter
	  
	  ClusterTempsAvg[i+j*Run.NCls] /= elapsed;

	  if(DebugPrint == true && j==0) //console printing...
	    printf("T[%s] %8.2f, ",Clusters[i].Name,ClusterTempsAvg[i+j*Run.NCls]);
	  
	  fwrite(&ClusterTempsAvg[i+j*Run.NCls],sizeof(float),1,MDstats); //-> cls temp(avg) as 1 float
	  
	}
	ClusterTempsAvg[i+j*Run.NCls] = 0.0f;  //reset the average clstemp /for all clusters
      }//end of clusters loop
      
      //reset the avg E
      AvgE[j] = 0;
      AvgEK[j]= 0;
      
    }//end of images loop
    
    if(DebugPrint == true)
      printf("\n");
 

    //---------- force part --------------
    AvgTipPos.x /= elapsed;
    AvgTipPos.y /= elapsed;
    AvgTipPos.z /= elapsed;

    //print the shit
    fwrite(&time,sizeof(float),1,MDforce);         //-> time(averaged) as 1 float
    fwrite(&AvgTipPos,sizeof(float3),1,MDforce);   //-> tip position(averaged) as 1 float3

    //forces on the tip clusters only!
    int cls;
    for(j=0;j<Run.MDimages;j++){
       
       for(i=0;i<Run.TipCls;i++){
	 cls = Run.TipClsIdx[i];
	 ClusterForce[cls+j*Run.NCls] = Scalar_f3_f( ClusterForce[cls+j*Run.NCls], 1.0f/elapsed ); //divide the average
	   fwrite(&ClusterForce[cls+j*Run.NCls],sizeof(float3),1,MDforce); //->cluster force(avg) as 1 float3
       }
       
       for(i=0;i<Run.NCls;i++) // reset all cluster forces for this image
	 ClusterForce[i+j*Run.NCls] = make_float3( 0.0f, 0.0f, 0.0f ); //reset the average
     }
    
    //reset the average positions
    AvgTipPos = make_float3(0.0f, 0.0f, 0.0f); //reset the average tip-pos
    AvgTipHld = make_float3(0.0f, 0.0f, 0.0f); //reset the average holder-pos  

    elapsed = 0;
    elapsedt= 0.0f;

  }

    
  
  
  void WriteCoords(char *basename, float4 *Charges, float3 *Speeds)
  {
    FILE *fp;
    int i, img;
    char filename[100];
    
    sprintf(filename,"%s.xyz",basename);
    fp = fopen(filename,"w");

    //print a vmd xyz file where each frame is the configuration of an image!

    fprintf(fp, "%i\n\n", Run.NAtoms);
    
    for(img=0;img<Run.MDimages;img++){
      for(i=0;i<Run.NAtoms;i++){
	
	fprintf(fp, "%s %f %f %f  %f %f %f\n", Types[System[i].TypeID].Name, 
		Charges[i+img*Run.NAtoms].x, Charges[i+img*Run.NAtoms].y, Charges[i+img*Run.NAtoms].z,
		Speeds[i+img*Run.NAtoms].x, Speeds[i+img*Run.NAtoms].y, Speeds[i+img*Run.NAtoms].z);
	
      }
      fprintf(fp,"\n\n");
    }
    fclose(fp);
    
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
*/
