// -*- c++ -*-

#include <stdio.h>
//#include <stdlib.h>
//#include <strings.h>
#define true 1
#define false 0


//take nloops for one image and print them to a single binary
//file for Cinema 4D loader script

int StartsWith(char strIn[], char strRef[]);



int main (int argc, char* argv[])
{

  char filename[100], nameout[100],aname[5];
  int i,j,k,nloops,img, natoms, nimages;
  int frames;
  FILE *fin, *fout,*ftmp;
  short3 pos;
  float3 speed;

  if(argc<4){
    printf("Insufficient parameters! Specify the base name, the amount of loops and the image number\n");
    return 0;
  }
  
  
  nloops = atoi(argv[2]);
  img  = atoi(argv[3]);
  frames = 0;

  sprintf(nameout , "%s.trajectory.img%i.ctj",argv[1],img); //open the output file
  

  //open the first file and get NATM and NIMG
  //all files should have the same amount of atoms, images and frames

  sprintf(filename,"%s.loop%i.ctj",argv[1],1);
  printf("Opening %s\n",filename);
  fin = fopen(filename,"rb");
  fout= fopen(nameout ,"wb");
  fread(&natoms, sizeof(int),1,fin);
  fread(&nimages,sizeof(int),1,fin);

  printf("NUMBER OF ATOMS: %i\n",natoms);
  
  //this counts the frames
  while(!feof(fin)){ //loop until EOF
    frames++;
    for(i=0;i<natoms;i++){  //loop on all atoms

      if((int)fread(aname,sizeof(char),5,fin) == 0)   //read the first atom
	break;

      fseek(fin,18*(img-1), SEEK_CUR);    //skip (img-1)(3short+3float)
      fread(&pos,sizeof(short3),1,fin);   //->position as float3
      fread(&speed,sizeof(float3),1,fin); //->speed as float3
      //printf("--%s-- %i %i %i -- %f %f %f\n",aname,pos.x,pos.y,pos.z,speed.x,speed.y,speed.z);
  
      fseek(fin,18*(nimages-img), SEEK_CUR);  //now skip (NIMG-img)(3short+3float)
    }
  }
  printf("TOTAL FRAMES IN FILE: %i\n",frames);
  fclose(fin);
  
  //ftmp = fopen("tmp.xyz","w");
  //fprintf(ftmp,"1\n\n");


  i=frames*nloops;
  fwrite(&natoms, sizeof(int),1, fout);
  fwrite(&i, sizeof(int),1, fout);
  
  for(j=0;j<natoms;j++){   //loop over all atoms
    for(i=1;i<=nloops;i++){   //loop over the files

      sprintf(filename,"%s.loop%i.ctj",argv[1],i); //loop filename
      fin = fopen(filename,"rb");  //open it
      fseek(fin,8,SEEK_CUR);  //skip the first 2 int      
      fseek(fin,(5+18*nimages)*j ,SEEK_CUR); //jump to the right atom (j)

      for(k=0;k<frames;k++){   //loop over all frames in the current file
	fread(aname,sizeof(char),5,fin);  //read name
	fseek(fin,18*(img-1), SEEK_CUR);  //skip non wanted images
	//read data for the wanted image
	fread(&pos,sizeof(short3),1,fin);   //->position as float3
	fread(&speed,sizeof(float3),1,fin); //->speed as float3
	if(j==0){
	  //printf("--%s-- %i %i %i -- %f %f %f\n",aname,pos.x,pos.y,pos.z,speed.x,speed.y,speed.z);
	  //fprintf(ftmp,"%s %f %f %f\n\n\n",aname,pos.x/100.0,pos.y/100.0,pos.z/100.0);
	}
	fseek(fin,18*(nimages-img), SEEK_CUR);  //skip the subsequent images
	fseek(fin,(5+18*nimages)*(natoms-1) ,SEEK_CUR); //skip natoms atoms/images

	if(k==0 && i==1) //only once per atom
	  fwrite(aname,sizeof(char),2,fout);
	fwrite(&pos,sizeof(short3),1,fout);

      }

      fclose(fin);
    }

  }

  fclose(fout);

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



