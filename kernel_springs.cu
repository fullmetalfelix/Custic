// -*- c++ -*-



//computes the coulomb forces between a and b
__device__ void SpringForce(float4 a, float4 b, float k, float Force[])
{
  //float3 Force;
  float r[3];
  
  
  r[0] = b.x - a.x;
  r[1] = b.y - a.y;
  r[2] = b.z - a.z;

  float distsq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  float dist = rsqrtf(distsq);
  float qq = a.w*b.w / distsq;

  #pragma unroll
  for(int c=0;c<3;c++)
    {
      r[c] *= dist;
      Force[c] -= r[c]*qq ;
    } 
  
}



//this kernel is called for the N cores
__global__ void Calc_Spring_Forces( float4 *pos, float3 *forces, float *Springs )
{
  
  __shared__ float4 cpos[BLOCKSIZE];
  __shared__ float4 spos[BLOCKSIZE];
  
  int idx = blockIdx.x*BLOCKSIZE + threadIdx.x;
  float myforce[3], r[3];
  short isok = 1;
  

  if(idx >= NAtoms_d)
    {
      idx = 0;
      isok = 0;
    }

  

  cpos[threadIdx.x] = pos[idx];
  
  spos[threadIdx.x] = pos[idx+NAtoms_d];
  
  
  myforce[0] = 0.0f;  myforce[1] = 0.0f;  myforce[2] = 0.0f;  //set the forces to 0
  
  //core->shell vector
  r[0] = cpos[threadIdx.x].x - spos[threadIdx.x].x;
  r[1] = cpos[threadIdx.x].y - spos[threadIdx.x].y;
  r[2] = cpos[threadIdx.x].z - spos[threadIdx.x].z;
  
  //float dist = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  //dist = sqrtf(dist);
  
  #pragma unroll
  for(int c=0; c<3; c++)
    {
      myforce[c] = r[c]*Springs[idx];//*SpringForceFactor_d;
    }

  //copy back from shared to global
  if(isok == 1)
    {
      forces[idx].x -= myforce[0];
      forces[idx].y -= myforce[1];
      forces[idx].z -= myforce[2];
      
      forces[idx+NAtoms_d].x += myforce[0];
      forces[idx+NAtoms_d].y += myforce[1];
      forces[idx+NAtoms_d].z += myforce[2];     

    }

  
}




