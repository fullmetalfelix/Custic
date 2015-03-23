// -*- c++ -*-


extern "C"
{
  



  

  //get the forces from the GPU positions
  //the result are copied back, processed //and copied on the gpu again!
  //forMD = true : no springs and shell forces goes on the core
  //forMD = false: shell forces stays on shell and the polarization spring is added (to the shell)
  int GetForces(int forMD, int DebugPrint)
  {
    
    dim3 dimBlock_os(BLOCKSIZE);
    dim3 dimGrid_os(Run.NAtoms / BLOCKSIZE + 1);
    //unsigned int t1;
  
    //cutCreateTimer(&t1);
    //cutStartTimer(t1);
    Calc_Forces<<<dimGrid_os, dimBlock_os>>>(Charges_d, AtomType_d, forces_d);
    cudaThreadSynchronize();
 
    //if(CopyBack == true)
    cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);  //copy back

      
    
    //cudaMemcpy(forces_d, forces_h, f3_NObj, cudaMemcpyHostToDevice); //send the processed forces on the gpu
     
    /*cutStopTimer(t1);

    printf("Forcer grid is %i\n",Run.NAtoms / BLOCKSIZE + 1);
    printf("Forces computed in %f ms\n",cutGetTimerValue(t1));
    FILE *fp = fopen("force1.out","w");
    int i;
    for(i=0;i<Run.NAtoms;i++)
      {
	fprintf(fp,"%8.5f %8.5f %8.5f - %8.5f \n",Charges_h[i].x,Charges_h[i].y,Charges_h[i].z,Charges_h[i].w);
	fprintf(fp,"%8.5f %8.5f %8.5f \n",forces_h[i].x,forces_h[i].y,forces_h[i].z);
      }
    fclose(fp);
    

    //cutDeleteTimer(t1);*/

    return true;
    
  }



  /*
  int GetForces_GPU_quiet()
  {

    unsigned int t1;
    int i,ii;
    dim3 dimBlock_os(BLOCKSIZE);
    dim3 dimGrid_os(Run.NAtoms*2 / BLOCKSIZE + 1);
    dim3 dimGrid_as(Run.NAtoms / BLOCKSIZE + 1);

    cutCreateTimer(&t1);
    cutStartTimer(t1);
    Calc_Forces<<<dimGrid_os, dimBlock_os>>>(Charges_d, AtomType_d, Springs_d, ShellParameters_d, forces_d);
    cudaThreadSynchronize();

    Sort_Forces<<<dimGrid_as, dimBlock_os>>>( Springs_d, forces_d );
    cudaThreadSynchronize();
 
    cudaMemcpy(forces_h, forces_d, f3_NObj, cudaMemcpyDeviceToHost);  //copy back
    cutStopTimer(t1);

    printf("Forces computed in %f ms\n",cutGetTimerValue(t1));
    

    FILE *fp = fopen("force2.out","w");
    for(i=0;i<Run.NAtoms;i++)
      {
	ii=i+Run.NAtoms;
	fprintf(fp,"%8.5f %8.5f %8.5f - %8.5f %8.5f %8.5f\n",forces_h[i].x,forces_h[i].y,forces_h[i].z,forces_h[ii].x,forces_h[ii].y,forces_h[ii].z);
      }
    fclose(fp);

    cutDeleteTimer(t1);
    return true;
  }
  */
  
  
}
