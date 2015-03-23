#include <stdio.h>
#include "defos.h"







void GetForce_onI(int i)
{
  
	int j,c;
	double dist[3], r;
  double force;

  for(j=0;j<Run.NAtoms;j++)  //loop on all atoms
    {
      
      if(i != j )
	{
	  //compute the distance
	  r = 0.0;
	  for(c=0;c<3;c++)
	    {
	      dist[c] = System[j].poscore[c] - System[i].poscore[c];
	      r += dist[c] * dist[c];
	    }
	  //r = sqrt(r);
	  force = Types[System[j].TypeID].ChargeIon * Types[System[i].TypeID].ChargeIon / r;
	  for(c=0;c<3;c++)
	    {
	      dist[c] /= r;
	      System[i].Force[c] = dist[c] * force;
	    }
	  //--------------------
	  
	  
	  
	}
      
      
      
      
    }
  
  
  
  
}

