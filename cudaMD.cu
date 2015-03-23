// -*- c++ -*-





//MD functions
extern "C"
{

  // Kb = 1.3806503e-23 SI
  // 1 uma = 1.660538782e-27 Kg
  // 1 eV  = 1.60217646e-19 J
  // 1 eV/(Åuma) = 0.964853382e18  m/s² = 0.964853382e-2 Å/fs² (acceleration)
  // 1 N = 0.602214179e7 uma Å/fs²
  // 1 Å/fs = 

  //positions are in Å
  //speeds are in Å/fs
  //forces are in eV/Å
  //time is in fs
  //masses are in uma

  void GetSystemTemp( void );


  void LeapFrog ( )
  {
    
    //this kernel seems faster with block 64!
    //SpdUpdate_plain<<<SlowGrid, SlowBlock>>>(forces_d, Speeds_d, Holders_d);  //update the velocities
    SpdUpdate<<<SlowGrid, SlowBlock>>>(forces_d, Speeds_d, Holders_d, Temperature_GPU);
    cudaThreadSynchronize();
    //*** now Speeds are the velocities at i+1/2; the thermostat is applied automatically
    
    GetSystemTemp(); //compute the temperature of the system (all images) and updates the value of gamma   

    //now we update the positions
    PosUpdate<<<SlowGrid, SlowBlock>>>(Charges_d, Speeds_d, Holders_d);  //update the positions
    cudaThreadSynchronize();
    //and the leapfrog is done!

  }


  

  
  /* Randomizes the speeds according to a set temperature, for all images of the system
     Remove the total average momentum from them so that Ptot = 0
     Then copy them to the gpu!
     NOTE: you are supposed to have the right speeds in the CPU memory!!!
   */
  void RandomizeSpeed( int Image )
  {
    int i,j;
    float x1,x2,y1,y2;
    float kt;
    //float T, T0;

    float3 Ptot;
    int jstart = 0, jend = Run.MDimages;
    if(Image>0){
      jstart = Image;
      jend = Image+1;
    }
    

    printf("Ramdomizing velocities (%i)...\n",Image);
    for(j=jstart; j<jend; j++) {
      //reset the Ptot
	Ptot = make_float3(0.0f,0.0f,0.0f);
      
      for(i=0;i<Run.NAtoms;i++) {

	x1 = (float)rand()/RAND_MAX;
	x2 = (float)rand()/RAND_MAX;
	//------------------------------ CAREFUL - READ T FROM INPUT! -------------
	kt = (float)sqrt( 2.0f*1.3806503e4 *  300.0f / (1.66053886f * Mass_h[i])) ;
	
	y1 = (float)(sqrt(-log(x1)) * cos( 2*3.1415f * x2 ) * kt  );
	y2 = (float)(sqrt(-log(x1)) * sin( 2*3.1415f * x2 ) * kt  );
	
	//printf("ktval = %f %f %f \n", kt,x1,x2);
	
	Speeds_h[i+j*Run.NAtoms].x = y1 * 1.0e-5;
	Speeds_h[i+j*Run.NAtoms].y = y2 * 1.0e-5;
	
	x1 = (float)rand()/RAND_MAX;
	x2 = (float)rand()/RAND_MAX;
	
	y2 = (float)(sqrt(-log(x1)) * sin( 2*3.1415f * x2 ) * kt  );
	Speeds_h[i+j*Run.NAtoms].z = y2 * 1.0e-5;

	Ptot.x += Speeds_h[i+j*Run.NAtoms].x * Mass_h[i];
	Ptot.y += Speeds_h[i+j*Run.NAtoms].y * Mass_h[i];
	Ptot.z += Speeds_h[i+j*Run.NAtoms].z * Mass_h[i];
      }
 
      //get the Ptot average
      Ptot.x /= 3.0f*Run.NAtoms;
      Ptot.y /= 3.0f*Run.NAtoms;
      Ptot.z /= 3.0f*Run.NAtoms;

      //remove the total average momentum from each atom
      for(i=0;i<Run.NAtoms;i++){
	Speeds_h[i+j*Run.NAtoms].x -= Ptot.x / Mass_h[i];
	Speeds_h[i+j*Run.NAtoms].y -= Ptot.y / Mass_h[i];
	Speeds_h[i+j*Run.NAtoms].z -= Ptot.z / Mass_h[i];
	//printf("%f %f %f\n",Speeds_h[i+j*Run.NAtoms].x,Speeds_h[i+j*Run.NAtoms].y,Speeds_h[i+j*Run.NAtoms].z);
      }

    }
    
    cudaMemcpy(Speeds_d, Speeds_h, f4_NObj, cudaMemcpyHostToDevice); //copy to gpu
    cutilCheckMsg("rndspd memcopy failed");
  }






}


 
