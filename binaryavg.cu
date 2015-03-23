// -*- c++ -*-

#include <stdio.h>
//#include <stdlib.h>
//#include <strings.h>
#define true 1
#define false 0
#define BUFFER 7000


//extract a given loop of a given image
int StartsWith(char strIn[], char strRef[]);

//usage: custicavg startloop endloops images fc

int main (int argc, char* argv[])
{

  char filename[100];
  int i,j,k, nimages;
  int iniloop, maxloop;
  FILE *fin, *fout, *fdissis;
  double dissi, *dissiimg;
  int ncls,line=0,comp;
  float time,ztip,avgforce[BUFFER], *loopforce;
  float tip[BUFFER]; //all tip positions
  float force[3];
  double avgdissisum = 0.0;

  if(argc<5){
    printf("Insufficient parameters!\n");
    return 0;
  }
  
  fdissis = fopen("dissis.txt","w");

  iniloop = atoi(argv[1]); printf("Loop start: %i\n",iniloop);
  maxloop = atoi(argv[2]); printf("Loop end  : %i\n",maxloop);
  nimages = atoi(argv[3]); printf("Number of images: %i\n",nimages);
  comp = atoi(argv[4]); printf("Force component: %i\n",comp);

  loopforce = (float*)malloc(sizeof(float)*nimages*BUFFER);
  dissiimg = (double*)malloc(sizeof(double)*nimages);
 

  for(j=0;j<BUFFER;j++){ //reset the average loop force
    avgforce[j] = 0.0f;
    tip[j] = 0.0f;
  }

  for(i=iniloop;i<=maxloop;i++){ //loop over the loops
    
    sprintf(filename,"force.MD.loop%i.out",i);
    printf("Opening force file loop #%i (%s)...\n",i,filename);
    
    fin = fopen(filename,"rb");
    if(fin == NULL){
      printf("Failed! File not found!\n");
      break;
    }
    
    dissi = 0.0;
    for(j=0;j<nimages;j++)
      dissiimg[j] = 0.0; //dissipation for this loop, for each image

    line = 0;

    //now loop until the end of file
    while (!feof(fin)){ //this read all records in the force curve
      
      fread(&time, sizeof(float),1,fin); //read the time
      if(feof(fin))
	break;
      

      fread(&time,sizeof(float),1,fin); //TIP X
      fread(&time,sizeof(float),1,fin); //TIP Y
      fread(&ztip,sizeof(float),1,fin); //TIP Z

      tip[line] = ztip;

      for(j=0;j<nimages;j++){ //loop over the images
	
	fread(&(force[0]),sizeof(float),1,fin); //read Fx
	fread(&(force[1]),sizeof(float),1,fin); //read Fy
	fread(&(force[2]),sizeof(float),1,fin); //read Fz
	fread(&time,sizeof(float),1,fin); //read tip  Temp
	fread(&time,sizeof(float),1,fin); //read surf Temp
	//printf("lread %i - %f!\n",line,force[comp]);
	
	//now fz is the force on the tip for IMAGE J
	avgforce[line] += force[comp]/nimages/(maxloop-iniloop+1); //sum force for all images
	loopforce[j*BUFFER + line] = force[comp];
	//printf("loopforce[%i]\n",j*BUFFER+line);
      }//end images loop

      
      line++;
    }//now all lines in the loop file were read
    printf("total lines: %i\n",line);

    for(k=0;k<nimages;k++){
      for(j=1;j<line;j++){ //compute the dissipation for THIS! loop, all images
	dissiimg[k] += (double)(loopforce[k*BUFFER+j]+loopforce[k*BUFFER+j-1])*
	  (tip[j]-tip[j-1])*0.5;
      }
      //printf("ASD %i\n",k);
      //dissiimg[k] /= (maxloop-iniloop+1);
      fprintf(fdissis,"%f ",dissiimg[k]);
    }

    for(k=0;k<nimages;k++)
      avgdissisum += dissiimg[k]/nimages;

    fprintf(fdissis,"\n"); //new line
    
    fclose(fin); //close the file and go to the beginning to open a new loop file
  } //END OF FILE LOOP!

  fout = fopen("avgloop.txt","w");
  for(j=0;j<line;j++){
    fprintf(fout,"%f %f\n",tip[j],avgforce[j]);
  }
  fclose(fout);
  
  dissi=0;
  for(j=1;j<line;j++){
    dissi += (double)(avgforce[j]+avgforce[j-1])*(tip[j]-tip[j-1])*0.5;
  }
  printf("avg dissipation: %f\n",dissi);
  
  //at this point avgdissisum is the sum of all individual dissi/nimages...
  avgdissisum /= (maxloop-iniloop+1); //divide by n loops
  printf("avg dissipation sum: %f\n",avgdissisum);

  free(dissiimg);
  free(loopforce);



  fclose(fdissis);
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


