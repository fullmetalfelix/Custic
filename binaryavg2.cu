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

  int i,j, nimages;
  FILE *fin;
  float time;
  float force[3];
 
  float3 tforce, avgforce,force1;

  if(argc<3){
    printf("Insufficient parameters!\n");
    return 0;
  }

  nimages = atoi(argv[2]);
  
  //printf("opening stats file: %s (%i images)\n",argv[1],nimages);
  fin = fopen(argv[1], "rb");

  // force file format is:
  // |time|xtip|ytip|ztip| |Fx1|Fy1|Fz1|T1t|T1s| |Fx2|Fy2|Fz2|T2t|T2s| ... |FxN|FyN|FzN|TNt|TNs|
  // |f32 |f32 |f32 |f32 | |f32|f32|f32|f32|f32| |f32|f32|f32|f32|f32| ... |f32|f32|f32|f32|f32|
  avgforce = make_float3(0.0f,0.0f,0.0f);
  i=0;

  while(!feof(fin)){
    
    //read the time
    fread(&time, sizeof(float),1,fin);
    
    //read the tip position
    fread(force, sizeof(float),3,fin);

    tforce = make_float3(0.0f,0.0f,0.0f);

    for(int j=0;j<nimages;j++){
      
      fread(&force1, sizeof(float3),1,fin); //read the force;
      tforce.x += force1.x;
      tforce.y += force1.y;
      tforce.z += force1.z;

      fread(force, sizeof(float),2,fin); //read the temperatures
    }
    
    avgforce.x += tforce.x/nimages;
    avgforce.y += tforce.y/nimages;
    avgforce.z += tforce.z/nimages;
    i++;
  }
  
  fclose(fin);

  printf("%f %f %f \n",avgforce.x/i, avgforce.y/i, avgforce.z/i);

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


