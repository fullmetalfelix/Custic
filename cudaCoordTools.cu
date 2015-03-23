// -*- c++ -*-


extern "C"
{
  void RandomizeSpeed( int ImageID );
  void ResetImage( int ImageID );
  void FirstHolder( void );


  /*
    put the tip in the position given in the input file 'position'
    referred to the sticker position 'tipstick' - all images
  */
  void TipReposition( int Move )
  {
    float3 dr;
 
    if(Move == true){
      dr.x = Run.TipXY[0]-Run.TipXY0[0];
      dr.y = Run.TipXY[1]-Run.TipXY0[1];    
      dr.z = Run.TipXY[2]-Run.TipXY0[2];

      printf("Moving the tip by: %f %f %f\n",dr.x,dr.y,dr.z);
      
      //move the tip and copy the positions to CPU
      MoveTip( dr );
      cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost);
      cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost); //holders too
      cudaMemcpy(HolderI_d, Holders_d, f4_NObj, cudaMemcpyDeviceToDevice); //and save the starting holders
    }

    TipPos.x = Run.TipXY[0];
    TipPos.y = Run.TipXY[1];
    TipPos.z = Run.TipXY[2];
    printf("The tip is now in: %f %f %f\n",TipPos.x,TipPos.y,TipPos.z);

    FirstHolder();
    TipHldDist.x = -TipPos.x + TipHld.x;
    TipHldDist.y = -TipPos.y + TipHld.y;
    TipHldDist.z = -TipPos.z + TipHld.z;
    printf("Tip-TipHolder distance: %f %f %f\n",TipHldDist.x,TipHldDist.y,TipHldDist.z);

  }

  /*
    starting from the input file configuration, put each image in a different point in space.
    the tip positions are taken from a vector that holdes the grid points, starting from offset.
   */
  void Grid_TipReposition(float3 *grid, int offset, int gridsize)
  {
    
    float3 dr[Run.MDimages];
    
    for(int i=0;i<Run.MDimages;i++)
      dr[i] = make_float3(0.0f,0.0f,0.0f);
    
    //store the displacements in the forces vector
    for(int i=offset; i<offset+Run.MDimages; i++){
      
      if(i==gridsize) //exit if we ran out of gridpoints to compute
	break;
      
      dr[i-offset].x = grid[i].x - Run.TipXY0[0];
      dr[i-offset].y = grid[i].y - Run.TipXY0[1];
      dr[i-offset].z = grid[i].z - Run.TipXY0[2];
      
    }
    
    
  }


  int IsTipAtom(int Index)
  {

    int c;

    for(c=0;c<Run.TipCls;c++)
      if(System[Index].ClusterID == Run.TipClsIdx[c])
	return true;

    return false;
  }


  //check if the tip came up with a chain attached.
  void CheckChain( int loop )
  {

    int i,j;
    float radius;
    float ztipi = 999.0f;
    int tipidx=Run.NAtoms-1;


    FILE *fp = fopen("changes.log","a");
    for(j=0;j<Run.MDimages;j++){

      for(i=0;i<Run.NAtoms;i++){  //find the lower tip atom
	if(Charges_h[i+j*Run.NAtoms].z <= ztipi && IsTipAtom(i)==true){
	  ztipi = Charges_h[i+j*Run.NAtoms].z;
	  tipidx = i;
	}
      }
      
      for(i=0;i<Run.NAtoms;i++){  //loop over the non tip atoms
	
	if(IsTipAtom(i)==false){
	  radius = pow(Charges_h[i+j*Run.NAtoms].x-Charges_h[tipidx+j*Run.NAtoms].x,2);
	  radius+= pow(Charges_h[i+j*Run.NAtoms].y-Charges_h[tipidx+j*Run.NAtoms].y,2);
	  radius+= pow(Charges_h[i+j*Run.NAtoms].z-Charges_h[tipidx+j*Run.NAtoms].z,2);
	  radius = sqrt(radius);
	  if(radius <= 3.0f){
	    printf("IMAGE %i has a chain or tip change (atom %i)\n",j,i);
	    fprintf(fp,"chain in loop %i image %i\n",loop,j); //prints a one in the cumulative chain log
	    ResetImage(j);
	    break;
	  }
	}
      }
      
    }//end of loop over images
	
    fclose(fp);

  }


  
  //reset an image to the starting coordinates (from chargeO) and randomize speeds
  void ResetImage(int ImageID)
  {
    int i;
    
    for(i=0;i<Run.NAtoms;i++){ 
      //take the coords from the restart written at the beginning of the loop
      Charges_h[i+ImageID*Run.NAtoms].x = ChargeO_h[i+ImageID*Run.NAtoms].x;
      Charges_h[i+ImageID*Run.NAtoms].y = ChargeO_h[i+ImageID*Run.NAtoms].y;
      Charges_h[i+ImageID*Run.NAtoms].z = ChargeO_h[i+ImageID*Run.NAtoms].z;
      //Speeds_h[i+ ImageID*Run.NAtoms] = make_float4(0.0f, 0.0f, 0.0f, Mass_h[i]);
    }
    cudaMemcpy(Charges_d,Charges_h,f4_NObj,cudaMemcpyHostToDevice); //copy to gpu
    RandomizeSpeed( ImageID );  //now reinitialize the speeds

    printf(" image %i has been reset.\n",ImageID);

  }





  //******************************************************************************************
  /*  GPU kernel for resetting the tip holder after a loop is completed.
      Holders for springed/fixed atoms are reloaded from HREF,
      completely frozen atoms will be reset too.
   */
  __global__ void ResetHolder_GPU( float4 *X, float4 *H, float4 *HREF, short *isTip)
  {

    float4 myX, myH, myREF;
    
    short4 asd; // x=imageindex   y=blockindex(in image)   z=atomindex   w=istip...
    asd.x = (blockIdx.x/FperIMG_d);
    asd.y = blockIdx.x - asd.x*FperIMG_d;
    asd.z = asd.y * BLOCKFAST + threadIdx.x;
    asd.x *= NAtoms_d;  // x=imageindex*Natoms
 
    short isok = (asd.z < NAtoms_d);
    asd.z *= isok;

    asd.w = isTip[asd.z];
    myX   = X[asd.z+asd.x];
    myH   = H[asd.z+asd.x];
    myREF = HREF[asd.z+asd.x];

    asd.w = ( (asd.w == 1) && (myH.w != 0) );
    if(asd.w == 1){ //if the atom is in the tip and is somehow constrained
      myH = myREF; //restore the original holder
    }
    asd.w = ( (asd.w == 1) && (myH.w < 0) );
    if(asd.w == 1){ //if the atom is fixed for real...
      myX.x = myREF.x;//put the atom on the holder
      myX.y = myREF.y;
      myX.z = myREF.z;
    }

    //if the atom exists, save it
    if(isok == 1){
      H[asd.z+asd.x] = myH;
      X[asd.z+asd.x] = myX;
    }
    
  }


  //reset the fixed part of the tip (the holder cluster) on the GPU
  void ResetHolder( )
  {

    ResetHolder_GPU<<<FastGrid, FastBlock>>>(Charges_d, Holders_d, HolderI_d, TipAtoms_d);
    cudaThreadSynchronize();

    //copy positions and holders back
    cudaMemcpy(Charges_h, Charges_d, f4_NObj, cudaMemcpyDeviceToHost);
    cudaMemcpy(Holders_h, Holders_d, f4_NObj, cudaMemcpyDeviceToHost);

  }
  
}
