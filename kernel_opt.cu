// -*- c++ -*-


//displace the atoms along the direction of the forces
__global__ void OptStep( float4 *pos, float3 *forces )
{
  
  int idx = threadIdx.x + blockIdx.x*BLOCKSIZE;

  if(idx >= NAtoms_d)
    return;

  float4 mypos = pos[idx];
  float3 myfor = forces[idx];

  float displ = myfor.x*myfor.x + myfor.y*myfor.y + myfor.z*myfor.z;
  displ = sqrtf(displ);
  
  if(displ > MaxDispl_d)
    displ = MaxDispl_d;

  displ = 1.0f/displ;

  mypos.x += myfor.x*displ;
  mypos.y += myfor.y*displ;
  mypos.z += myfor.z*displ;
  
   pos[idx] = mypos;
   
}