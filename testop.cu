// -*- c++ -*-





//test functions functions
//extern "C"
//{


  float3 *vel_d;
  float *mas_d;
  float *res_d;
  unsigned int *clsid_d;

  void test_alloc()
  {

    
    
    cudaMalloc((void**)&vel_d, f3_NAtoms);
    cudaMalloc((void**)&mas_d, f1_NAtoms);
    cudaMalloc((void**)&res_d, /*sizeof(float)*(Run.NAtoms/BLOCKSIZE+1)*/f1_NAtoms);
    cudaMalloc((void**)&clsid_d, u1_NAtoms);


    unsigned int cid[Run.NAtoms];
    float res[(Run.NAtoms)];
    int i;

    for(i=0;i<Run.NAtoms;i++)
      {
	cid[i] = (unsigned int)System[i].ClusterID;
	res[i] = 0.0f;
      }
    
    cudaMemcpy(mas_d, Mass_h, f1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(res_d, res, f1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(clsid_d, cid, u1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(vel_d, Speeds_h, f3_NAtoms, cudaMemcpyHostToDevice);


  }

  

  void test_free()
  {

    cudaFree(vel_d);
    cudaFree(mas_d);
    cudaFree(res_d);
    cudaFree(clsid_d);

  }

  __device__ unsigned int count = 0;
  __global__ void SumVels( float3 *speeds, float *masses, unsigned int *cID, volatile float *results)
  {
    __shared__ bool isLastBlockDone;
    __shared__ float3 svel[BLOCKSIZE];
    __shared__ float  smas[BLOCKSIZE];
    __shared__ unsigned int scid[BLOCKSIZE];
    __shared__ float psum[BLOCKSIZE];

    int idx = threadIdx.x + blockIdx.x * BLOCKSIZE;  //atom index
    int tdx = threadIdx.x; //thread index and cluster
    int isok = 1;
    unsigned int value;

    if(idx >= NAtoms_d)
      {
	isok = 0;
	idx = 0;
      }

    svel[tdx] = speeds[idx];
    smas[tdx] = masses[idx];
    scid[tdx] = cID[idx];
    psum[tdx] = 0.0f;

    __syncthreads();
    
    
    for(int i=0;i<BLOCKSIZE;i++) //sum the speeds
      if( scid[i] == tdx )
      	psum[tdx] += isok * 0.5f * smas[i] * ( svel[i].x*svel[i].x + svel[i].y*svel[i].y + svel[i].z*svel[i].z ) * 103.64269f;
    
    
    __syncthreads();

    
    if(isok == 1)
      results[idx] = (idx==tdx)*(psum[tdx]/0);
    /*
    __syncthreads();

    
    __threadfence();

    if(tdx == 0) //thread 0 sends the increment
      {

	//__threadfence();

	value = atomicInc(&count, gridDim.x);

	isLastBlockDone = (value == (gridDim.x - 1));


      }

    //__threadfence();
    __syncthreads(); //all wait for thread 0 who sent the increment and was fenced
    
    
    if (isLastBlockDone) //if this is the last block
      {
	psum[tdx] = 0.0f; //each thread resets...
	for(int i=0;i<gridDim.x;i++)
	  if(i*BLOCKSIZE+tdx < NAtoms_d)
	    psum[tdx] += results[i*BLOCKSIZE+tdx];//each thread loads a partial result and sum

	__syncthreads();

	results[tdx] = psum[tdx];

	if(tdx == 0)
	  {
	    count = 0;
	    results[0] = 0;
	  }
	 
      }
    
    */
  }


 
  void testtemp()
  {
    
    dim3 blk(BLOCKSIZE);
    dim3 grd(Run.NAtoms / BLOCKSIZE + 1);

    unsigned int t1;
    float res[Run.NAtoms];
    
    test_alloc();


    cutCreateTimer(&t1);
    cutStartTimer(t1);
    SumVels<<< grd, blk >>>(vel_d, mas_d, clsid_d, res_d);
    cudaThreadSynchronize();

    cudaMemcpy(res, res_d, f1_NAtoms, cudaMemcpyDeviceToHost);
    cutStopTimer(t1);

    
    int i,j;
    float temps[BLOCKSIZE];

    for(i=0;i<BLOCKSIZE;i++)
      temps[i] = 0.0f;
    
    
    for(i=0;i<grd.x;i++)
      for(j=0;j<BLOCKSIZE;j++)
	if(i*BLOCKSIZE + j < Run.NAtoms)
	  temps[j] += res[i*BLOCKSIZE + j];
    
    /*
    for(i=0;i<5+0*Run.NAtoms;i++)
      printf("res[%i] = %f \n",i,Speeds_h[i].x );
    */
   
    /*
    for(j=0;j<3*Run.NCls;j++)
      printf("T[%i] = %f - %f - %f\n",j,res[j],ClusterKins[j],ClusterTemps[j]);
    */

    printf("test computed in %f ms\n",cutGetTimerValue(t1));
    printf("count is %i\n",count);

    cutDeleteTimer(t1);
    
    //cudaMemcpy(&count,&i,sizeof(unsigned int),cudaMemcpyHostToDevice);

    test_free();

  }


 void testtemp_cpu()
  {
    
    unsigned int t1;
 
    cutCreateTimer(&t1);
    cutStartTimer(t1);
    
    int i,j,c,x,y;
    float force[3];
    float r[3], dist, qq;
    float4 p;

    for(i=0;i<2*Run.NAtoms;i++)
      {
	force[0]  = 0.0f; force[1]  = 0.0f;force[2]  = 0.0f;
	for(j=0;j<2*Run.NAtoms;j++)
	  {
	    if(i!=j)
	      {
		r[0] = Charges_h[j].x - Charges_h[i].x;
		r[1] = Charges_h[j].y - Charges_h[i].y;
		r[2] = Charges_h[j].z - Charges_h[i].z;

		dist = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
		dist = sqrtf(dist);
		
		if( (j!=i+Run.NAtoms) && (j!=i-Run.NAtoms) )
		  {
		    qq = dist*dist*dist;
		    qq = Charges_h[i].w * Charges_h[j].w *14.3996368420f / qq;
                    #pragma unroll
		    for(int c=0;c<3;c++)
		      force[c] -= r[c]*qq;
		  }
		
		if( (j >= Run.NAtoms) && (i>=Run.NAtoms)  )
		  {
		    x = AtomType_h[i];
		    y = AtomType_h[j];
		    p =  ShellParameters_h[ max(x,y)*(max(x,y)+1)/2 + min(x,y) ];
		    qq = (-p.x*expf(-dist/p.y)/p.y + p.z*p.w*powf(dist,p.w-1)) / dist;
                    #pragma unroll
		    for(int c=0;c<3;c++)
		      force[c] += r[c]*qq ;
		  }

	      }
	    
	    
	    
	  }
      }


    
    cutStopTimer(t1);
    
    
    printf("test computed in %f ms\n",cutGetTimerValue(t1));
    cutDeleteTimer(t1);
    
    //cudaMemcpy(&count,&i,sizeof(unsigned int),cudaMemcpyHostToDevice);
    

  }




//}
