extern "C"
{
  //utilities and wrappers
  float3 Sum_f3_f3(float3 a, float3 b)
  {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
  }
  float4 SumFloat4_3(float4 a, float4 b)
  {
    if(a.w != b.w)
      printf("WARNING! SumFloat4_3 summing vectors with different W.\n");
    return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, a.w);
  } 

  float4 SumFloat4_Float3(float4 a, float3 b)
  {
    return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, a.w);
  } 

  float3 Dif_f4_f4(float4 a, float4 b)
  {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
  }
  float3 Dif_f3_f3(float3 a, float3 b)
  {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
  }
  float Scalar_f3(float3 a, float3 b)
  {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }
  float Scalar_f4(float3 a, float3 b)
  {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }
  float3 Scalar_f3_f(float3 a, float c)
  {
    return make_float3(a.x*c, a.y*c, a.z*c);
  }
}