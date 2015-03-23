// -*- c++ -*-


extern "C"
{

  void AllocateAtoms();
  int  AllocatePairPot( PairPot *CPUPairs );
  int  AllocateTemperature( void );
  void AllocateSymbols();
  void AllocateEnergy();
  int  AllocateAvgForces();
  void SetThermoBounds( );

  //define standard array sizes
  void Alloc_Sizes()
  {
    f1_NObj   =     Run.NAtoms * sizeof(float)  * Run.MDimages;
    f4_NObj   =     Run.NAtoms * sizeof(float4) * Run.MDimages; //for position/charge, speed/mass, holder/spring
    f3_NObj   =     Run.NAtoms * sizeof(float3) * Run.MDimages; //for forces 

    f1_NAtoms =     Run.NAtoms * sizeof(float); //for masses
    f3_NAtoms =     Run.NAtoms * sizeof(float3);
    f4_NAtoms =     Run.NAtoms * sizeof(float4);

    i3_NAtoms =     Run.NAtoms * sizeof(int3); //for ion fixers

    s1_NAtoms =     Run.NAtoms * sizeof(short);
    s3_NAtoms =     Run.NAtoms * sizeof(short3) * Run.MDimages;

    f3_ShellPars =  MAX_PPOT * sizeof(float3); //size of pairpot parameters table

    printf("Common sizes initialized.\n");

    //define the grid/block sizes
    Run.BperIMG = ceil((float)Run.NAtoms / BLOCKSIZE);
    printf("Blocks per image: %i\n",Run.BperIMG);

    Run.FperIMG = ceil((float)Run.NAtoms / BLOCKFAST);
    printf("Blocks per image: %i (fast kernels)\n",Run.FperIMG);

    FastGrid = max((Run.FperIMG * Run.MDimages),1);
    FastBlock= BLOCKFAST;
    SlowGrid = max((Run.BperIMG * Run.MDimages),1);
    SlowBlock= BLOCKSIZE;    
    dim3 gridt(Run.MDimages,2);
    TGrid = gridt;

  }




//CPURun has the operation parameters, System has the atoms and Types the type descriptors
  int AllocateGPU( RunParams CPURun,  Atom *CPUSystem, AtomType *CPUTypes, 
		   ClusterDef *CPUClusters, PairPot *CPUPairs )
  {
    //int i,j;
    cudaError_t err;
  
    //copy the starting config from the cpu version
    Run = CPURun; 
    System = CPUSystem;
    Types = CPUTypes;
    Clusters = CPUClusters;

    //choose a device
    err = cudaSetDevice(Run.DeviceID);
    cutilCheckMsg("alloc setdev failed");
    if(err != 0 )
      {
	printf("cudaSetDevice ERROR %i \n", err);
	return false;
      }
  
    //define standard array sizes
    Alloc_Sizes();
    cutilCheckMsg("alloc sizes failed");
    
    //allocate forces vector
    forces_h = (float3*) malloc(f3_NObj);
    cudaMalloc((void**)&forces_d, f3_NObj);

    printf("Force vectors allocated!\n");cutilCheckMsg("f?");
    //**************************************

    AllocatePairPot( CPUPairs );cutilCheckMsg("pairs?");
    
    AllocateTemperature();cutilCheckMsg("temp?");
    AllocateAvgForces();cutilCheckMsg("avgf?");
    AllocateAtoms();cutilCheckMsg("atoms?");
    SetThermoBounds( );cutilCheckMsg("thermob?");
    //AllocateEnergy();
    
    printf("  Total number of degrees of freedom: %i\n", NDegFree);
 
 
    //copy constant memory static symbols
    AllocateSymbols();
    

    //allocate resources for Conj grad linear system solver
    //AllocateConjLin();
    
   
    
    printf("Data copied to the GPU.\n");
    cutilCheckMsg("really?");
    return true;

  }

 
  
  //allocate vectors with atomic positions
  void AllocateAtoms()
  {

    int i,j;
    
    //allocate the positions - for each image
    Charges_h = (float4*) malloc(f4_NObj);
    ChargeO_h = (float4*) malloc(f4_NObj);
    ChargeT_h = (float4*) malloc(f4_NObj);
    ChargeINIT = (float4*)malloc(f4_NObj);
   
    cudaMalloc((void**)&Charges_d, f4_NObj);
    cudaMalloc((void**)&ChargeT_d, f4_NObj);

    //holders and springs for fixed atoms
    Holders_h = (float4*) malloc(f4_NObj);   //one holder for each atom - All images
    cudaMalloc((void**)&Holders_d, f4_NObj); //not needed but maybe improves mem access
    cudaMalloc((void**)&HolderI_d, f4_NObj);
    cudaMalloc((void**)&HolderR_d, f4_NAtoms);
    

    //masses and inverse masses
    Mass_h = (float*) malloc(f1_NAtoms);
    iMass_h = (float*) malloc(f1_NAtoms);
    cudaMalloc((void**)&iMass_d, f1_NAtoms);
    cudaMalloc((void**)&Mass_d, f1_NAtoms);
  
    //types
    AtomType_h = (short*) malloc(s1_NAtoms);
    cudaMalloc((void**)&AtomType_d, s1_NAtoms);
    
    //Speeds
    Speeds_h = (float4*) malloc(f4_NObj);
    cudaMalloc((void**)&Speeds_d, f4_NObj);
 
    //tip atoms
    TipAtoms_h = (short*)malloc(s1_NAtoms);
    cudaMalloc((void**)&TipAtoms_d, s1_NAtoms);

    
    Run.NTipAtoms = 0;
    //--- create atomic properties vectors ---
    for(i=0;i<Run.NAtoms;i++){	
      Mass_h[i] = Types[System[i].TypeID].Mass;  //mass
      iMass_h[i] = 1.0f/Mass_h[i];               //inverse mass
      AtomType_h[i] = (short)System[i].TypeID;   //type index
      
      TipAtoms_h[i] = 0;
      for(j=0;j<Run.TipCls;j++){  //check if the atom is in the tip
	if( System[i].ClusterID == Run.TipClsIdx[j] ){
	  Run.NTipAtoms++;
	  TipAtoms_h[i] = 1;
	  break;
	}
      }
      
    }
    cudaMemcpy( Mass_d,  Mass_h,  f1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy( iMass_d, iMass_h, f1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy( AtomType_d, AtomType_h, s1_NAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy( TipAtoms_d, TipAtoms_h, s1_NAtoms, cudaMemcpyHostToDevice);
    //--------------------------------------------------------

    //create positions and speeds
    for(j=0;j<Run.MDimages;j++)
      for(i=0;i<Run.NAtoms;i++)	{
	Speeds_h[i + j*Run.NAtoms] = make_float4(0.0f, 0.0f, 0.0f,Mass_h[i]);	
	Charges_h[i + j*Run.NAtoms] = make_float4( System[i].poscore[0], System[i].poscore[1], System[i].poscore[2], 
						   Types[System[i].TypeID].Charge);
	Holders_h[i + j*Run.NAtoms] = make_float4( System[i].poscore[0], System[i].poscore[1], System[i].poscore[2],
						   Clusters[System[i].ClusterID].spring);
      }
    cudaMemcpy( Speeds_d,  Speeds_h,  f4_NObj, cudaMemcpyHostToDevice); //vx,vy,vz,m
    cudaMemcpy( Charges_d, Charges_h, f4_NObj, cudaMemcpyHostToDevice); //x,y,z,q
    cudaMemcpy( Holders_d, Holders_h, f4_NObj, cudaMemcpyHostToDevice); //cx,cy,cz,spring
    //------------------------------------------------------
    
    //relative positions for HOLDERS
    float4 HolderRel[Run.NAtoms];
    for(i=0;i<Run.NAtoms;i++){
      HolderRel[i].x = Holders_h[i].x - Holders_h[0].x;
      HolderRel[i].y = Holders_h[i].y - Holders_h[0].y;
      HolderRel[i].z = Holders_h[i].z - Holders_h[0].z;
    }
    HolderRel[0] = make_float4(0.0f,0.0f,0.0f,0.0f);
    cudaMemcpy( HolderR_d, HolderRel, f4_NAtoms, cudaMemcpyHostToDevice);

    
    printf("Atoms and charges allocated!\n");

  }

  int AllocateTemperature()
  {
    //allocate GPU temperature vectors - use bigger blocks as low shared is needed
    Temperature_CPU = (float*)malloc(2*Run.MDimages*sizeof(float));
    cudaMalloc((void**)&Temperature_GPU,2*Run.MDimages*sizeof(float));
 
    for(int i=0;i<2*Run.MDimages;i++){
      Temperature_CPU[i] = 0.0f;
    }
    cudaMemcpy(Temperature_GPU, Temperature_CPU, 2*Run.MDimages*sizeof(float), cudaMemcpyHostToDevice);
    
    printf("Temperature vectors allocated.\n");
    return true;
  }

  int AllocateAvgForces()
  {
    //allocate GPU temperature vectors - use bigger blocks as low shared is needed
    TipForce_CPU = (float3*)malloc(Run.MDimages*sizeof(float3));
    cudaMalloc((void**)&TipForce_GPU,Run.MDimages*sizeof(float3));
    
    for(int i=0;i<Run.MDimages;i++)
      TipForce_CPU[i] = make_float3(0.0f,0.0f,0.0f);
    cudaMemcpy(TipForce_GPU, TipForce_CPU, Run.MDimages*sizeof(float3), cudaMemcpyHostToDevice);
    
    printf("Average Forces vectors allocated.\n");
    return true;
  }

  int AllocatePairPot( PairPot *CPUPairs )
  {
    
    //allocate the shell parameter table CONSTANT MEM
    PairPot p;
    int i,j;

    ShellParameters_h = (float3*)malloc(f3_ShellPars);
    //cudaMalloc((void**)&ShellParameters_d, f4_ShellPars);

    if(f3_ShellPars > sizeof(float3)*MAX_PPOT)
      {
	printf("FATAL! Too many pair potentials! (AllocatePairPot)\n");
	return false;
      }

    for(i=0;i<Run.NTypes;i++)
      for(j=0;j<Run.NTypes;j++){
	p = CPUPairs[i*Run.NTypes + j];
	ShellParameters_h[i*Run.NTypes + j] = make_float3(p.Value[0], 1.0f/p.Value[1], p.Value[2]);
	ShellParameters_h[j*Run.NTypes + i] = make_float3(p.Value[0], 1.0f/p.Value[1], p.Value[2]);
      }
    
    cudaMemcpyToSymbol( ShellParameters_d, ShellParameters_h, f3_ShellPars, 0, cudaMemcpyHostToDevice);

    for(i=0;i<Run.NTypes;i++)
      for(j=0;j<Run.NTypes;j++)
	printf("pairpot %i %i -> %f %f %f\n",i,j,ShellParameters_h[i*Run.NTypes + j].x,ShellParameters_h[i*Run.NTypes + j].y,ShellParameters_h[i*Run.NTypes + j].z);

    printf("Pair potentials allocated in CONSTANT memory.\n");
    
    return true;

  }


  //allocate symbols in constant memory
  void AllocateSymbols ( )
  {

    //copy constant memory static symbols
    cudaMemcpyToSymbol( NImages_d, &Run.MDimages, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( BperIMG_d, &Run.BperIMG, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( FperIMG_d, &Run.FperIMG, sizeof(int), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol( NAtoms_d, &Run.NAtoms, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( NTipAtoms_d, &Run.NTipAtoms, sizeof(unsigned short), 0, cudaMemcpyHostToDevice);
    
    cudaMemcpyToSymbol( MDstep_d, &Run.MD_Step, sizeof(float), 0, cudaMemcpyHostToDevice);

    float q = Run.MD_Step;
    cudaMemcpyToSymbol( MDQstep_d, &q, sizeof(float), 0, cudaMemcpyHostToDevice);
    

    cudaMemcpyToSymbol( TipBound_d, &Run.TipBound, sizeof(short), 0, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&TipPosZ_d, sizeof(float));

  }

  //finds the start/end points of tip & surface thermostats
  void SetThermoBounds( )
  {
    int isin;

    //tip thermo
    for(int i=0;i<Run.NAtoms;i++){
      isin = 0;
      for(int c=0;c<Run.ThermoTipCls;c++)
	if(System[i].ClusterID == Run.ThermoTip[c]){
	  ThermoBounds_h.x = i;
	  isin = 1;
	  break;
	}
      if(isin == 1)
	break;
    }
    for(int i=ThermoBounds_h.x; i<Run.NAtoms;i++){
      isin = 0;
      for(int c=0;c<Run.ThermoTipCls;c++)
	if(System[i].ClusterID == Run.ThermoTip[c])
	  isin++; //isin is 0 if the atom is not in any of the thermo tip clusters
      if(isin == 0){
	ThermoBounds_h.y = i; //index of the first atom in the non thermo cluster
	break;
      }
      
    }

    //surf thermo
    for(int i=0;i<Run.NAtoms;i++){
      isin = 0;
      for(int c=0;c<Run.ThermoSurfCls;c++)
	if(System[i].ClusterID == Run.ThermoSurf[c]){
	  ThermoBounds_h.z = i;
	  isin = 1;
	  break;
	}
      if(isin == 1){
	//printf("AAA! first thermosurf found %i c%i(%s)\n",i,System[i].ClusterID,Clusters[System[i].ClusterID].Name);
	break;
      }
    }
    
    ThermoBounds_h.w = Run.NAtoms; //upper limit in case the cluster is at the end of the list
    for(int i=ThermoBounds_h.z; i<Run.NAtoms;i++){
      isin = 0;
      for(int c=0;c<Run.ThermoSurfCls;c++)
	if(System[i].ClusterID == Run.ThermoSurf[c]){
	  isin++; //isin is 0 only if itz not in the thermosurface
	}    
      if( (isin == 0) || (i == Run.NAtoms-1) ){ //if we are out of the thermo region or at the end of the atom list...
	ThermoBounds_h.w = i; //index of the first atom in the non thermo cluster
	break;
      }
    }

    printf("Tip  thermostat on atoms %i-%i\n",ThermoBounds_h.x,ThermoBounds_h.y);
    printf("Surf thermostat on atoms %i-%i\n",ThermoBounds_h.z,ThermoBounds_h.w);

    float cpl[2], tset[2];
    cpl[0] = 1.0f/Run.TCoupTip;
    tset[0]= Run.TempTip;
    cpl[1] = 1.0f/Run.TCoupSurf;
    tset[1]= Run.TempSurf;
    
    cudaMemcpyToSymbol( ThermoBounds_d, &ThermoBounds_h, sizeof(short4), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( Tset_d, tset, 2*sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( Tcup_d, cpl , 2*sizeof(float), 0, cudaMemcpyHostToDevice);

  }



  //******************************************************************************************************
  //*** DEALLOCATOR ***
  int FreeGPU()
  {
    
    printf("deallocating...\n");
    //free the atomic positions
    free(Charges_h); cudaFree(Charges_d);
    free(ChargeT_h); cudaFree(ChargeT_d);
    free(ChargeO_h); free(ChargeINIT);

    free(Holders_h); cudaFree(Holders_d);
    cudaFree(HolderI_d); cudaFree(HolderR_d);

    //free the speeds
    free(Speeds_h); cudaFree(Speeds_d);
 
    //free the atomic properties
    free(forces_h); cudaFree(forces_d);
    
    free(Mass_h); cudaFree(Mass_d);
    free(iMass_h);cudaFree(iMass_d);

    free(AtomType_h); cudaFree(AtomType_d);
    free(TipAtoms_h); cudaFree(TipAtoms_d);

    free(ShellParameters_h); //cudaFree(ShellParameters_d);


    //free the temperature stuff
    free(Temperature_CPU); cudaFree(Temperature_GPU);
 
    //free the avg tip forces
    free(TipForce_CPU); cudaFree(TipForce_GPU);
    
    cudaFree(TipPosZ_d);

    System = 0;

    //free(AvgE); free(AvgEK); free(Etot);

    
    printf("GPU resources were deallocated.\n");
  
    return true;
  }





}
