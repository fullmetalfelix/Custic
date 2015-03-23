// -*- c++ -*-

#include <string.h>

extern "C"
{

  //void DoOptStep_Fragment(float step, int DebugPrint);   //this is in cudainterface_CPU
  

  float GetKinetic(int j)
  {
    int i;
    float result = 0.0f;
    //float hv2;

    //take it from the clusterkins
    for(i=0;i<Run.NCls;i++)
      result += ClusterKins[i+j*Run.NCls];

    return result;


  }




}
