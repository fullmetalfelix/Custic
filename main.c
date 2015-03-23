
#include <stdio.h>
#include "defos.h"
//#include <cudef.ch>

//extern int AllocateGPU( RunParams CPURun, Atom *CPUSystem, AtomType *CPUTypes );
//extern int FreeGPU();


int main (int argc, char *argv[])
{

  //first thing, check the input file
  if(CheckInputFile(argc, argv) == false )
  {
    printf("FATAL! Bad input file!\n");
    return 15;
  }

  if(ReadInput(argv[1]) == false)
    {
      printf("stop\n");
      return false;
    }


  printf("Allocating GPU resources...\n");
  AllocateGPU( Run, System, Types, Clusters, PairPots);
  
  //xyzMake("traj.xyz");
  

  //Steepest(true);
  //GetEnergy(true, true);
  //GetForces(true,true);
  //LineSearch(true);
  //Steepest_LS(true);
  
  //DoOptStep(true);
  //BFGS(true);
  //ConjGradOptimize(true),
  
  MainProcess();

  printf("Main process done!\n");

  //testmatmul();
  //testmatvec();

  DeAllocate(); printf("CPU dealloc done!\n");
  FreeGPU();

  printf("%i\n", sizeof(float));
  printf("%i\n", sizeof(int));
  printf("%i\n", sizeof(short));

  return true;
}




//Deallocate all the resources
int DeAllocate()
{

  //free(AtomSystem);

  if(Types != 0)
    free(Types);
  if(Clusters != 0)
    free(Clusters);
  if(System != 0)
    free(System);

}



int Map(int i, int j)
{ return i*Run.NTypes+j;}

