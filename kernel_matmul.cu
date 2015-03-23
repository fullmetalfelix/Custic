// -*- c++ -*-



// perform A*B -> C
__global__ void MatMatMul(float *A, float *B, float *C)
{

  //SQUARE MATRICES!!!
  __shared__ float tileA[MBLOCK*MBLOCK];
  __shared__ float tileB[MBLOCK*MBLOCK];
  
  float result = 0.0f;
  int i,j;
 
  int myIdxj = blockIdx.x*MBLOCK + threadIdx.x;
  int myIdxi = blockIdx.y*MBLOCK + threadIdx.y;

  int uj = MBLOCK * blockIdx.x;      //skip bdim*bidx columns
  int ui = MBLOCK * blockIdx.y * nd;  //skip bdim*bidx lines  
  
  for(i=0; i< nd / MBLOCK; i++)
    {
      
      //each thread loads one guy in shared
      tileA[threadIdx.y*MBLOCK + threadIdx.x] = A[ui + threadIdx.x + threadIdx.y*nd + MBLOCK*i];
      tileB[threadIdx.y*MBLOCK + threadIdx.x] = B[uj + threadIdx.x + threadIdx.y*nd + MBLOCK*i*nd];
      __syncthreads();
      
      for(j=0; j<MBLOCK; j++)
	{
	  result += tileA[threadIdx.y*MBLOCK + j] * tileB[j*MBLOCK + threadIdx.x];
	}

      __syncthreads();
      

    }

    C[myIdxi*nd + myIdxj] = result;  
    
}

//multiply A and v and store the result in w - RUN FOR N THREADS?
__global__ void MatVecMul(float *A, float *v, float *w)
{


  __shared__ float tileA[BLOCKSIZE];
  __shared__ float tilev[BLOCKSIZE];
  
  float result = 0.0f;
  //float final = 0.0f;
  int i,j;
 
  int myIdxi = blockIdx.x;

  for(i=0; i< LinSize_d / BLOCKSIZE; i++)
    {

      //each thread loads one guy in shared
      tileA[threadIdx.x] = A[blockIdx.x*LinSize_d + threadIdx.x +  BLOCKSIZE*i];
      tilev[threadIdx.x] = v[threadIdx.x + BLOCKSIZE*i];
      __syncthreads();
      
      for(j=0; j<BLOCKSIZE; j++)
	  result += tileA[j] * tilev[j];

      __syncthreads();
      
    }

  if(threadIdx.x == 0)
    w[myIdxi] = result;  

}

//do b - k*Ap  -- k is alpha in the const mem
__global__ void GetResidue(float *b, float *A, float *p, float* r)
{


  __shared__ float tileA[BLOCKSIZE];
  __shared__ float tilep[BLOCKSIZE];
  
  float result = 0.0f;
  //float final = 0.0f;
  int i,j;
 
  int myIdxi = blockIdx.x;

  for(i=0; i< LinSize_d / BLOCKSIZE; i++)
    {
      
      //each thread loads one guy in shared
      tileA[threadIdx.x] = A[blockIdx.x*LinSize_d + threadIdx.x +  BLOCKSIZE*i];
      tilep[threadIdx.x] = p[threadIdx.x + BLOCKSIZE*i];
      __syncthreads();
      
      for(j=0; j<BLOCKSIZE; j++)
	  result += tileA[j] * tilep[j];

      __syncthreads();
      
    }

  if(threadIdx.x == 0)
    r[myIdxi] = b[myIdxi] - result;  

}


/*
__global__ void ScalProd(float *v, float *w, float *r)
{
  
  __shared__ float tilev[BLOCKSIZE];
  __shared__ float tilew[BLOCKSIZE];
  
  float result = 0.0f;
  int i,j;


  for(i=0; i< LinSize_d / BLOCKSIZE; i++)
    {
      tilev[threadIdx.x] = v[threadIdx.x + BLOCKSIZE*i];
      tilew[threadIdx.x] = w[threadIdx.x + BLOCKSIZE*i];
      __syncthreads();

      for(j=0; j<BLOCKSIZE; j++)
	result += tilev[j] * tilew[j];

      __syncthreads();
      
    }
  
  if(threadIdx.x == 0)
    r[blockIdx.x] = result;
  

}
*/
