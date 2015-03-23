#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "defos.h"

extern int ProbeDevice( int );



int ReadInput(char *filename)
{
  FILE *fin;


  fin = fopen(filename, "r");


  //read the atomic types
  if(ReadTypes(fin) == false)
    return false;
  
  //read the interaction parameters
  if(ReadPairs(fin) == false)
    return false;
  
  //read the clusters definitions
  if(ReadClusters(fin) == false)
    return false;

  //read the atomic coordinates
  if(ReadAtoms(fin) == false)
    return false;
  
  
  if(ReadMD(fin) == false)  //read the MD box if present
    return false;
  if(ReadTipDyn(fin) == false)
    return false;
  
  //sort the atoms; tip clusters first
  int i,j,c,tc, istip;
  printf("Sorting atoms...\n");
  j = 0;
  //first put the tip thermostat clusters
  for(c=0;c<Run.ThermoTipCls;c++){
    tc = Run.ThermoTip[c];
    for(i=0;i<Run.NAtoms;i++){ //loop over all atoms to find the one in the tip
      if(tSystem[i].ClusterID == tc){//if the atom is in the tip cluster, put it in System
	System[j] = tSystem[i];
	j++;
      }
    }
  }//tip thermo clusters written
  //printf("tipthermo done, j is %i\n",j);

  //then the other tip clusters
  for(c=0;c<Run.TipCls;c++){ //loop over the tip clusters
    tc = Run.TipClsIdx[c];
    istip = 0;
    for(i=0;i<Run.ThermoTipCls;i++)
      if(tc == Run.ThermoTip[i])
	istip = 1; //istip goes to 1 if this cluster is found in the thermotip
    if(istip == 1)
      continue; //and if it was, skip it
    
    //if it wasnt in the thermotip, write its atoms
    for(i=0;i<Run.NAtoms;i++){ //loop over all atoms to find the one in the tip
      if(tSystem[i].ClusterID == tc){//if the atom is in the tip cluster, put it in System
	System[j] = tSystem[i];
	j++;
      }
    }
  }//FULL TIP WRITTEN BY NOW
  //printf("fulltip done, j is %i\n",j);
  Run.TipBound = j;

  //write the surface thermostat clusters
  for(c=0;c<Run.ThermoSurfCls;c++){
    tc = Run.ThermoSurf[c];
    for(i=0;i<Run.NAtoms;i++){ //loop over all atoms to find the one in the tip
      if(tSystem[i].ClusterID == tc){//if the atom is in the tip cluster, put it in System
	System[j] = tSystem[i];
	j++;
      }
    }
  }//surf thermo clusters written
  //printf("thermosurf done, j is %i\n",j);


  //then the other surface clusters
  for(c=0;c<Run.NCls;c++){ //loop over the tip clusters
    istip=-1; //assume itz not tip
    for(tc=0;tc<Run.TipCls;tc++){ //loop over the tip clusters
      if(c == Run.TipClsIdx[tc]){  //if it is tip, skip it
	istip = 1;
	break;
      }
    }
    for(i=0;i<Run.ThermoSurfCls;i++) //now check if it is in the thermosurf
      if(c == Run.ThermoSurf[i])
	istip = 1; //istip goes to 1 if this cluster is found in the thermotip
    if(istip == 1)
      continue; //and if it was, skip it
    
    //if it wasnt in the tip nor thermosurf, write its atoms
    for(i=0;i<Run.NAtoms;i++){ //loop over all atoms to find the one in the tip
      if(tSystem[i].ClusterID == c){//if the atom is in the tip cluster, put it in System
	System[j] = tSystem[i];
	j++;
      }
    }
  }//FULL TIP WRITTEN BY NOW
  printf("j is %i\n",j);
  
  //all done
  free(tSystem);
  
  fclose(fin);

  printf("Input file processed correctly.\n");
  return true;

}



int ReadAtoms(FILE *fin)
{

  int i = 0, j, isNull, cIdx;
  int words;
  int typeOK = false;
  char line[READBUFFER], tmpstr[20];
  float core[3], shll[3];
  
  if( FindString(fin,"<Atoms>") == false )
    {
      printf("FATAL! Missing Atoms block!\n");
      return false;
    }

  //count the atoms
  while(StartsWith(fgets(line, sizeof(line), fin), "<End>") != true)
    {
      isNull = false;
      for(j=0;j<Run.NCls;j++)
	if(StartsWith(line, Clusters[j].Name) == true)
	  {
	    isNull = true;
	    break;
	  }
      if(isNull == false){
	
	if(StartsWith(line,"#")==true) //avoid empty lines and comments
	  continue;
	words = CountWords(line);
	if(words >= 4)	
	  i++;
      }
    }
  if(i == 0)
    {
      printf("FATAL! No atom was found in the atoms block!\n");
      return false;
    }

  Run.NAtoms = i; printf("  Number of atoms: %u\n", Run.NAtoms);
  sSize = Run.NAtoms*( Run.NAtoms-1 )/2;
  System = (Atom*) malloc(Run.NAtoms * sizeof(Atom));
  tSystem = (Atom*) malloc(Run.NAtoms * sizeof(Atom));
  

  cIdx = 0; i = 0;
  fseek(fin,0,SEEK_SET); //rewind all
  FindString(fin,"<Atoms>");
  while(StartsWith(fgets(line, sizeof(line), fin), "<End>") != true)
    {
      isNull = false; 
      for(j=0;j<Run.NCls;j++)   //check if the line is actually a cluster name
	if(StartsWith(line, Clusters[j].Name) == true)
	  {
	    cIdx = j;
	    //printf("changing cluster to: %s\n",Clusters[cIdx].Name);
	    isNull = true;
	    break;
	  }

      if(isNull == false)
	{
	  //add the atom!
	  words = CountWords(line);
	  
	  if(StartsWith(line,"#")==true) // 
	    continue;

	  //case where we have name poscore
	  if(words >= 4){
	    sscanf(line, "%s %f %f %f",&tmpstr,&(tSystem[i].poscore[0]),&(tSystem[i].poscore[1]),
		   &(tSystem[i].poscore[2]));
	  } else  {
	    printf("FATAL! Insufficient parameters for atom #%i.\n",i);
	    return false;
	  }
	  //Associate the read name (tmpstr) to a type
          typeOK = false;
	  for(j=0;j<Run.NTypes;j++)
	    if(StartsWith(tmpstr,Types[j].Name) == true)
	      {
		tSystem[i].TypeID = j;
                typeOK = true;
		break;
	      }
	  if(typeOK == false)
	    {
	      printf("FATAL! The type of this (%s) atom was not defined!\n",tmpstr);
	    }
	  tSystem[i].ClusterID = cIdx;
	  	      
	  i++;
	}
    }

  fseek(fin,0,SEEK_SET); //rewind all
  return true;
}



int ReadClusters(FILE *fin)
{
  int i,j;
  char line[READBUFFER],tmpstr[20], tmpstr2[20], strval[READBUFFER];
  int words;
  
  
  if(FindString(fin,"<Clusters>") == false)
    {
      printf("FATAL! Missing Clusters definition box!\n");
      return false;
    }

  //the cluster defs was found
  //now count the clusters
  i = 0;
  while(StartsWith(fgets(line, sizeof(line), fin), "<End>") != true)
    i++;
  if(i == 0) {
    printf("FATAL! No cluster was defined in the clusters box!\n");
    return false;      
  }
  
  Run.NCls = i; printf("  Number of clusters: %i\n",Run.NCls);
  Clusters = (ClusterDef*) malloc(Run.NCls * sizeof(ClusterDef));
  fseek(fin,0,SEEK_SET); //rewind all
  FindString(fin,"<Clusters>");
  for(i=0;i<Run.NCls;i++) {
    fgets(line, sizeof(line), fin);
    words = CountWords(line);
    if(words != 3){
      printf("FATAL! Insufficient Cluster parameters in the <Cluster> box.\n");
      return false;
    }

    //case name, free|spring, 
    sscanf(line,"%s %s %f",&Clusters[i].Name, &strval, &Clusters[i].spring);
    if( StartsWith(strval, "free") == true){
      Clusters[i].hfix = 0;
      Clusters[i].spring = 0.0f;
    }
    else if( StartsWith(strval, "spring") == true)
      Clusters[i].hfix = 1;
    else if( StartsWith(strval, "fix") == true){
      Clusters[i].hfix = 1;
      Clusters[i].spring = -1.0f;  //fixed atoms have negative spring const
    }
    else{
      printf("ERROR! Invalid parameter in cluster definition.\n");
      return false;
    }
  }
  
  fseek(fin,0,SEEK_SET); //rewind all
  return true;

  
}


int ReadPairs(FILE *fin)
{

  int i,j,c;
  char line[READBUFFER],tmpstr[20], tmpstr2[20];
  int found,found2;
  PairPot *p;

  //for each pair of types there should be a <XXX-YYY> block

  //allocate the pair definition array
  //PairPots = (PairPot*)malloc(sizeof(PairPot)*Run.NTypes*(Run.NTypes-1)/2);
  PairPots = (PairPot*)malloc(sizeof(PairPot)*Run.NTypes*Run.NTypes);

  for(i=0;i<Run.NTypes;i++)
    for(j=0;j<=i;j++)
      {
	sprintf(tmpstr,"<%s-%s>",Types[i].Name,Types[j].Name);
	sprintf(tmpstr2,"<%s-%s>",Types[j].Name,Types[i].Name);
	//printf("CHECKING! %s \n",tmpstr);
	fseek(fin,0,SEEK_SET); //rewind all
	found = FindString(fin,tmpstr); found2 = true;
	if(found == false)
	  {
	    printf("INFO! Interaction box %s not found...\n",tmpstr);
	    fseek(fin,0,SEEK_SET); //rewind all
	    found2 = FindString(fin,tmpstr2);
	    if(found2 == false){
	      printf("INFO! Interaction box %s also not found...\n",tmpstr2);
	    }
	  }
	if(found2 == false)
	  {
	    printf("FATAL! Missing Interaction box <%s-%s> from the input file!\n",Types[i].Name,Types[j].Name);
	    return false;      
	  }
	//here we found the pair <ij>
	//now parse to the <End> and COUNT the interaction types/parameters
	//printf("Reading pair: %s\n",tmpstr);
	p = &PairPots[i*Run.NTypes + j];//&PairPots[Map(i,j)];
	while( StartsWith(fgets(line, sizeof(line), fin),"<End>") != true )
	  {
	    if(StartsWith(line,"buckingham") == true)
	      {
		p->Type = 1;
		sscanf(line,"%s %f %f %f %f",&tmpstr,&(p->Value[0]),&(p->Value[1]),&(p->Value[2]),&(p->Value[3]));
		//printf("pairpot read: %f %f %f %f\n",PairPots[i*Run.NTypes + j].Value[0],PairPots[i*Run.NTypes + j].Value[1],PairPots[i*Run.NTypes + j].Value[2],PairPots[i*Run.NTypes + j].Value[3]);
		//PairPots[j*Run.NTypes + i] = PairPots[i*Run.NTypes + j];
	      }
	    else
	      {
		printf("FATAL! Unrecognized pair potential label (%s).\n",tmpstr);
		return false;
	      }
	    
	  }

      }

  


  fseek(fin,0,SEEK_SET); //rewind all
  
  return true;
}



int ReadTypes(FILE *fin)
{
  char line[READBUFFER],tmpstr[READBUFFER];
  int i = 0;

  //check if the Types box exists
  if(FindString(fin,"<Types>") == false)
    {
      printf("FATAL! Missing Types box from the input file!\n");
      return false;      
    }

  //first read to the <End> and count how many types
  while(fgets(line, sizeof(line), fin) != NULL) {
    if(StartsWith(line,"<End>") == true)
      break;
    i++;
  }

  if(i == 0) {
    printf("FATAL! No atomic type is defined!\n");
    return false;
  }
  
  Run.NTypes = i;
  Types = (AtomType*)malloc(sizeof(AtomType)*Run.NTypes);   //allocate the types vector
  printf("  Number of atom types: %i\n",Run.NTypes);
  fseek(fin,0,SEEK_SET); //rewind all
  FindString(fin,"<Types>");
  for(i=0;i<Run.NTypes;i++){   //read $NTypes lines
    fgets(line, sizeof(line), fin);
    sscanf(line,"%s %f %f",&Types[i].Name, &Types[i].Mass, &Types[i].Charge );
  }
  fseek(fin,0,SEEK_SET); //rewind all

  return true;
}


//read the parameters for the simulation from the input file *fin
int ReadMD(FILE *fin)
{
  
  char line[READBUFFER],tmpstr[READBUFFER], strval[READBUFFER], strval2[READBUFFER];
  float tval, fval;
  int i,ival,words;

  //check if the parameters box exists
  if(FindString(fin,"<MD>") == false) {
    fseek(fin,0,SEEK_SET); //rewind all
    printf("ERROR! No Molecular Dynamics box was given!\n");
    return false;      
  }

  //if the box is there... the cursor is at the first line in the box

  Run.DeviceID = 0;
  if(ParseToEnd(fin, "gpuID", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %i",&tmpstr,&Run.DeviceID);  //split the line
  printf("  CUDA Device number: %i ",Run.DeviceID);
  if(ProbeDevice(Run.DeviceID) != true)
    return false;


  Run.MD_SaveXYZ = false;
  if(ParseToEnd(fin, "xyz", line,"<MD>") == true){   //parse until <End>
    sscanf(&line,"%s %s %i",&tmpstr,&strval, &Run.MD_XYZfreq);  //split the line
    if( (StartsWith(strval,"true") == true) || (StartsWith(strval,"TRUE") == true) ){
      Run.MD_SaveXYZ = true;
      printf("  MD xyz movie will be saved.\n");
    } 
  }else{
    printf("  MD xyz movie will NOT be saved.\n");
  }

  
  Run.MD_equiloop = false; //default= equilibrate only first???
  if(ParseToEnd(fin, "equiloop", line,"<MD>") == true){   //parse until <End>
    sscanf(&line,"%s %s",&tmpstr,&strval);  //split the line
    if( StartsWith(strval,"loops") == true || 
	StartsWith(strval,"all") == true ){
      Run.MD_equiloop = true;
      printf("  Equilibration will be performed before each loop.\n");
    }
    else if( StartsWith(strval,"once") == true ){
      Run.MD_equiloop = false;
      printf("  Equilibration will be performed in the beginning only.\n");
    }
  }else{
    printf("  Equilibration will be performed in the beginning only.\n");
  }


  Run.MD_reinit = false;
  if(ParseToEnd(fin, "reinit", line,"<MD>") == true){   //parse until <End>
    sscanf(&line,"%s %s",&tmpstr,&strval);  //split the line
    if( StartsWith(strval,"true") == true ){
      Run.MD_reinit = true;
      printf("  System will be reinitialized after every loop.\n");
    }
    else if( StartsWith(strval,"once") == true ){
      Run.MD_reinit = false;
      printf("  System will NOT be reinitialized after every loop.\n");
    }
  }else{
    printf("  System will NOT be reinitialized after every loop.\n");
  }

  //*** ***** RESTART FILE MANAGEMENT ***** ***
  Run.MD_reload = false; Run.MD_reload_num = 0; //default values: no reload, start from 0
  if(ParseToEnd(fin, "resume", line,"<MD>") == true){   //parse until <End>
    words = CountWords(line);
    if(words != 3){
      printf("ERROR! Insufficient parameters for resume directive!\n");
      return false;
    }
    sscanf(&line,"%s %s %i",&tmpstr,&strval,&Run.MD_reload_num);  //split the line
    if( (StartsWith(strval,"true") == true) || (StartsWith(strval,"TRUE") == true) ){
      Run.MD_reload = true;
      printf("  MD will resume from the restart point (loop number %i)!\n",Run.MD_reload_num);
    }
    else{
      printf("  MD will resume ?%s?\n",strval);
    }
  }
  //read the name of the restartfile
  if(Run.MD_reload == true){
    if(ParseToEnd(fin, "restartfile", line,"<MD>") == true){
      words = CountWords(line);
      if(words != 2){
	printf("ERROR! Restart file not specified correctly!\n");
	return false;
      }
      sscanf(&line,"%s %s",&tmpstr,&Run.MD_reloadfile);
      printf("  Restart file: %s\n", Run.MD_reloadfile);
    }
    else{
      printf("ERROR! Restart file not specified.\n");
      return false;
    }
  }
  //********************************************


  //*** TIP Thermostat ***
  if(ParseToEnd(fin, "thermotip", line,"<MD>") == true){   //parse until <End>
    words = CountWords(line);
    if(words < 4){
      printf("ERROR! Tip thermostat was not defined correctly.\n");
      return false;
    }
    TrimString(line);
    
    GetWord(line,strval,2); sscanf(strval,"%f",&Run.TempTip); //second word is temperature
    GetWord(line,strval,3); sscanf(strval,"%f",&Run.TCoupTip);//third word is coupl const
    
    //for each other word...
    for(i=4;i<=words;i++){
      GetWord(line,strval,i); //get the word (this should be the cluster name)
      ival = GetCluster(strval);
      if(ival<0)return false;
      
      Run.ThermoTip[i-4] = ival;
    }
    Run.ThermoTipCls = words-3;
    printf("  Tip  thermostat (%fK %f) on cluster: ",Run.TempTip,Run.TCoupTip);
    for(i=0;i<words-3;i++)
      printf("%s ",Clusters[Run.ThermoTip[i]].Name);
    printf("\n");
    
  }else{
    printf("ERROR! No thermostat for the tip was defined.\n");
    return false;
  }
  //******************************************************************************
  //*** SURF Thermostat ***
  if(ParseToEnd(fin, "thermosurf", line,"<MD>") == true){   //parse until <End>
    words = CountWords(line);
    if(words < 4){
      printf("ERROR! Surface thermostat was not defined correctly.\n");
      return false;
    }
    
    TrimString(line);
    
    GetWord(line,strval,2); sscanf(strval,"%f",&Run.TempSurf);
    GetWord(line,strval,3); sscanf(strval,"%f",&Run.TCoupSurf);    
    
    //for each other word...
    for(i=4;i<=words;i++){
      GetWord(line,strval,i); //get the word
      ival = GetCluster(strval);
      if(ival<0)return false;
      
      Run.ThermoSurf[i-4] = ival;
    }
    Run.ThermoSurfCls = words-3;
    printf("  Surf thermostat (%fK %f) on cluster: ",Run.TempSurf,Run.TCoupSurf);
    for(i=0;i<words-3;i++)
      printf("%s ",Clusters[Run.ThermoSurf[i]].Name);
    printf("\n");
  }else{
    printf("ERROR! No thermostat for the surface was defined.\n");
    return false;
  }    
  //******************************************************************************

  Run.MDimages = 1;
  if(ParseToEnd(fin, "images", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %i",&tmpstr, &Run.MDimages);  //split the line
  printf("  Number of system images: %i\n", Run.MDimages);

  Run.MD_Step = 1.0f;
  if(ParseToEnd(fin, "timestep", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %f",&tmpstr,&Run.MD_Step);  //split the line
  printf("  MD time step: %f ps\n", Run.MD_Step);
    
  Run.MD_NSteps = 200;
  if(ParseToEnd(fin, "steps", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %i",&tmpstr,&Run.MD_NSteps);  //split the line
  printf("  MD number of steps: %i\n", Run.MD_NSteps);

  Run.MD_ESteps = 5000;
  if(ParseToEnd(fin, "equisteps", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %i",&tmpstr,&Run.MD_ESteps);  //split the line
  printf("  MD number of steps (equilibration): %i\n", Run.MD_ESteps);


  Run.MD_Stats = 100;
  if(ParseToEnd(fin, "stats", line,"<MD>") == true)   //parse until <End>
    sscanf(&line,"%s %i",&tmpstr,&Run.MD_Stats);  //split the line
  printf("  MD statistics every: %i steps\n", Run.MD_Stats);
 
  fseek(fin,0,SEEK_SET); //rewind all

  return true;
}

int ReadTipDyn(FILE *fin)
{

  char line[READBUFFER],tmpstr[READBUFFER], strval[READBUFFER], strval2[READBUFFER];
  float tval, fval;
  int i,j;
  char *pch;

  float tipd[5];
  
  Run.TipDyn = 0;

  //check if the parameters box exists
  if(FindString(fin,"<Tip>") == false){
    fseek(fin,0,SEEK_SET); //rewind all
    printf("FATAL! No Tip description was given!\n");
    return false;      
  }
  
  //get the tip-clusters list
  if(ParseToEnd(fin, "tipcluster", line,"<Tip>") == true){   //parse until <End>
    
    i = CountWords(line);
    if(i<2) {
      printf("ERROR! Where do I apply the tip dynamics?\n");
      Run.TipDyn = 0;
      fseek(fin,0,SEEK_SET);
      return false;
    }
    i--;
    Run.TipCls = i; //total number of tip clusters
    
    printf ("  Tip clusters: ");
    pch = strtok (line," ");
    pch = strtok (NULL, " ");
    printf (" (holder %s) ",pch);
    for(j=0;j<Run.NCls;j++)
      if(StartsWith(pch,Clusters[j].Name) == true)
	Run.TipClsIdx[0] = j;
    pch = strtok (NULL, " ");

    for(i=0;i<Run.TipCls-1;i++) {
      printf (" %s",pch);

      for(j=0;j<Run.NCls;j++)
	if(StartsWith(pch,Clusters[j].Name) == true)
	  Run.TipClsIdx[i+1] = j;
      pch = strtok (NULL, " ");
    }
      
  } else {
    printf("ERROR! No valid tip description is present.\n");
    Run.TipDyn = 0;
    fseek(fin,0,SEEK_SET);
    return false;
  }
  //------------------------------------------------------------------------

  //read the tip-stickpoint
  if(ParseToEnd(fin, "tipstick", line,"<Tip>") == true)   //parse until <End>
    {
      sscanf(&line,"%s %f %f %f",&tmpstr,&Run.TipXY0[0],&Run.TipXY0[1],&Run.TipXY0[2]);
      printf("The tip stickpoint is at: %8.5f %8.5f %8.5f \n",Run.TipXY0[0],Run.TipXY0[1],Run.TipXY0[2]);
    }
  else
    {
      printf("FATAL! The tip stickpoint is not specified.\n");
      return false;
    }
  //------------------------------------------------------------------------
  //read the initial tip position
  if(ParseToEnd(fin, "position", line,"<Tip>") == true)   //parse until <End>
    {
      sscanf(&line,"%s %f %f %f",&tmpstr,&Run.TipXY[0],&Run.TipXY[1],&Run.TipXY[2]);
      printf("The tip will be repositioned at: %8.5f %8.5f %8.5f\n",Run.TipXY[0],Run.TipXY[1],Run.TipXY[2]);
    }
  else
    {
      printf("FATAL! The initial tip position is not specified.\n");
      return false;
    }
  //------------------------------------------------------------------------
   //read the resetter option
  Run.ChainReset = 0;
  if(ParseToEnd(fin, "chaincheck", line,"<Tip>") == true)   //parse until <End>
    {
      sscanf(&line,"%s %s",&tmpstr,&strval);
      if( StartsWith(strval,"true")==true ){
	printf("Geometry will be reset if chains or tip changes are present.\n");
	Run.ChainReset = 1;
      }else{
	printf("No action if chains or tip changes happens!\n");
      }
    }
  else
    {
      printf("FATAL! What should I do if a chain is formed?.\n");
      return false;
    }
  //------------------------------------------------------------------------ 



  //read the type of tip motion
  if(ParseToEnd(fin, "dynamics", line,"<Tip>") == true)   //parse until <End>
    {
      sscanf(&line,"%s %s",&tmpstr,&strval);

      if(StartsWith(strval, "static") == true){      //static
	printf("Tip dynamics: STATIC\n");
	Run.TipDyn = 0;
      }
      else if(StartsWith(strval, "linear") == true){ //linear motion
	printf("Tip dynamics: LINEAR\n");
	Run.TipDyn = 1;

	//read the LINEAR motion parameters

	//DIRECTION
	if(ParseToEnd(fin, "tip-direction", line,"<Tip>") == true){   //parse until <End>
	  sscanf(&line,"%s %f %f %f",&tmpstr,&Run.TipDynPars[0],&Run.TipDynPars[1],&Run.TipDynPars[2]);
	  tval = Run.TipDynPars[0]*Run.TipDynPars[0] + Run.TipDynPars[1]*Run.TipDynPars[1] +
	    Run.TipDynPars[2]*Run.TipDynPars[2];
	  tval = sqrt(tval);
	  for(i=0;i<3;i++)
	    Run.TipDynPars[i]/=tval;
	  printf("   direction: %8.5f %8.5f %8.5f\n",Run.TipDynPars[0],Run.TipDynPars[1],Run.TipDynPars[2]);
	}
	else{
	  printf("FATAL! Linear dynamics requires a direction!\n");
	  return false;
	}

	//SPEED
	if(ParseToEnd(fin, "tip-speed", line,"<Tip>") == true){   //parse until <End>
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[3]);
	  printf("       speed: %8.5f ang/ps\n",Run.TipDynPars[3]);
	}
	else{
	  printf("FATAL! Linear dynamics requires a speed!\n");
	  return false;
	}
	
      }//--------------------------------------------------------------
      else if(StartsWith(strval, "harmonic") == true){ //harmonic motion
	printf("Tip dynamics: HARMONIC\n");
	Run.TipDyn = 2;

	//read the HARMONIC motion parameters

	//FREQUENCY
	if(ParseToEnd(fin, "tip-frequency", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[0]);
	  printf("   frequency: %8.5f MHz\n",Run.TipDynPars[0]);
	}
	else{
	  printf("FATAL! Harmonic dynamics requires a frequency!\n");
	  return false;
	}	

	//AMPLITUDE
	if(ParseToEnd(fin, "tip-amplitude", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[1]);
	  printf("   amplitude: %8.5f ang\n",Run.TipDynPars[1]);
	}
	else{
	  printf("FATAL! Harmonic dynamics requires an amplitude!\n");
	  return false;
	}	

	//START POSITION
	if(ParseToEnd(fin, "tip-loopcut", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[2]);
	  printf("     loopcut: %8.5f ang\n",Run.TipDynPars[2]);
	}
	else{
	  Run.TipDynPars[2] = 0.0f;
	  printf("INFO! Harmonic dynamics starting time not found, using default (0.0).\n");
	  fseek(fin,0,SEEK_SET);
	  return true;
	}
	//AMOUNT OF LOOPS
	if(ParseToEnd(fin, "tip-nloops", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[3]);
	  printf("  # of loops: %i loops\n",(int)ceil(Run.TipDynPars[3]));
	}
	else{
	  Run.TipDynPars[3] = 1.0f;
	  printf("  # of loops: %i loop (default)\n",(int)ceil(Run.TipDynPars[3]));
	  fseek(fin,0,SEEK_SET);
	  return true;
	}
	

      }//--------------------------------------------------------------
      else if(StartsWith(strval, "TR") == true){ //harmonic motion
	printf("Tip dynamics: TORSIONAL \n");
	Run.TipDyn = 3;

	//read the HARMONIC motion parameters

	//FREQUENCY
	if(ParseToEnd(fin, "tip-frequency", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[0]);
	  printf("   frequency: %8.5f MHz\n",Run.TipDynPars[0]);
	}
	else{
	  printf("FATAL! Harmonic dynamics requires a frequency!\n");
	  return false;
	}	

	//AMPLITUDE
	if(ParseToEnd(fin, "tip-amplitude", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[1]);
	  printf("   amplitude: %8.5f ang\n",Run.TipDynPars[1]);
	}
	else{
	  printf("FATAL! Harmonic dynamics requires an amplitude!\n");
	  return false;
	}	

	//START POSITION
	if(ParseToEnd(fin, "tip-loopcut", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[2]);
	  printf("     loopcut: %8.5f ang\n",Run.TipDynPars[2]);
	}
	else{
	  Run.TipDynPars[2] = 0.0f;
	  printf("INFO! Harmonic dynamics starting time not found, using default (0.0).\n");
	  fseek(fin,0,SEEK_SET);
	  return true;
	}
	//AMOUNT OF LOOPS
	if(ParseToEnd(fin, "tip-nloops", line,"<Tip>") == true){   
	  sscanf(&line,"%s %f",&tmpstr,&Run.TipDynPars[3]);
	  printf("  # of loops: %i loops\n",(int)ceil(Run.TipDynPars[3]));
	}
	else{
	  Run.TipDynPars[3] = 1.0f;
	  printf("  # of loops: %i loop (default)\n",(int)ceil(Run.TipDynPars[3]));
	  fseek(fin,0,SEEK_SET);
	  return true;
	}
	
      }//--------------------------------------------------------------


      else{
	printf("ERROR! The tip dynamics specified is not valid.\n");
	fseek(fin,0,SEEK_SET);
	Run.TipDyn = 0;
	return true;
      }
    }
  





  for(i=0;i<5;i++)
    Run.TipDynPar2[i] = 0.0f;
  if(ParseToEnd(fin, "2dynamics", line,"<Tip>") == true)   //parse until <End>
    {
      i = CountWords(line);
      if(i < 6)
	{
	  printf("ERROR! The tip dynamics specified is not valid.\n");
	  fseek(fin,0,SEEK_SET);
	  return true;
	}
      sscanf(&line,"%s %s %f %f %f %f %f",&tmpstr,&strval,&tipd[0],&tipd[1],&tipd[2],&tipd[3],&tipd[4]);
      if(StartsWith(strval, "2harmonic") == true)
	{
	  //parameters are f (MHz) - Amp (nm) - offset (ang) - start point (us) - loops
	  printf("  2nd Harmonic tip dynamics will be added.\n");
	  for(i=0;i<5;i++)
	    Run.TipDynPar2[i] = tipd[i];
	  Run.TipDyn = 10;
	}

      else
	{
	  printf("ERROR! The tip dynamics specified is not valid.\n");
	  fseek(fin,0,SEEK_SET);
	  Run.TipDyn = 0;
	  return true;
	}
    }

  
  fseek(fin,0,SEEK_SET);
  return true;
}





int CheckInputFile(int argc, char *argv[])
{

  FILE *fin;

  if(argc<2)
    {
      printf("You need to specify an input file!\n");
      return false;
    }

  if( (fin=fopen(argv[1],"r")) == false)
    {
      printf("Invalid input file!\n");fclose(fin);
      return false;
    }
  //fclose(fin);
  printf("Input file exists!\n");
  return true;
}


int FindString(FILE *fin, const char strRef[] )
{
  char line[READBUFFER];

  while(fgets(line, sizeof(line), fin) != NULL )
      if(StartsWith(line,strRef) == true )
	return true;      
  return false;

}


int CountWords(char line[])
{
  
  int i,s=0,w=0;
  int l = strlen(line);
  int lastspace = true;

  //printf("len of line is %i characters \n",l);

  for(i=0;i<l;i++)
    {
      if((line[i] == ' ') || (line[i] == '\t'))
	{
	  if(lastspace != true)
	    {
	      s++;
	      lastspace = true;
	    }
	}
      else
	{
	  if(lastspace == true)
	    w++;
	  lastspace = false;
	  
	}
      
    }
  //printf("numbuh of spaces %i\n",s);
  //printf("numbuh of words  %i\n",w);

  return w;

}


/*
  Parse the lines until the first <End>
  return true if a line starting with refStr was found
  and set the outline pointer to the full line string
*/
int ParseToEnd(FILE *fin, char refStr[], char outline[], char Repo[])
{
  
  char line[READBUFFER];
  int found = false;

  do
    {
      fgets(line, sizeof(line), fin);           //read a line
      //printf("line was: %s�%s\n",line,refStr);
      
      if(StartsWith(line, refStr) == true)      //if itz the line we were looking for...
	{
	  //printf("FOUND THE PARSED STRING!\n");
	  memcpy(outline, line, READBUFFER);
	  found = true;
	  break;
	}

      if(StartsWith(line,"<End>") == true)	//check if the block is ended!
	{
	  //printf("FOUND THE END OF BLOCK!\n");
	  //fseek(fin,0,SEEK_SET);
	  break;
	}

    }while(!feof(fin));
  
  fseek(fin,0,SEEK_SET); //rewind all
  FindString(fin, Repo); //reposition the cursor in the beginning of the block!


  return found;
}

void GetWord(char line[],char out[], int word)
{

  if(word == 1){
    sscanf(line,"%s",out);
    return;
  }
  
  int i,j,s;
  s=1;
  for(j=0;j<strlen(line);j++){
    if(line[j] == ' ')
      s++;
    if(s==word)
      break;
  }
  
  sscanf(&line[j+1],"%s",out);
}



int GetCluster(char name[])
{
  int i;
  for(i=0;i<Run.NCls;i++){
    if(StartsWith(name,Clusters[i].Name) == true)
      break;
  }

  if(i>=Run.NCls){
    printf("ERROR! Cluster %s does not exist.\n",name);
    return -1;
  }

  return i;
}




void TrimString (char str[])
{
  char tmp[READBUFFER];
  int l = strlen(str);
  int lastspace = false;
  int i,j,start;

  //find the first non space character
  for(i=0;i<l;i++){
    if(str[i] != ' '){
      start = i;
      break;
    }
  }
  
  lastspace = false;
  j=0;
  for(i=start;i<l;i++){
    if((str[i] == ' ') && (lastspace == false)){
      lastspace = true;
      tmp[j] = str[i];
      j++;
    }
    if( str[i] != ' ' ){
      lastspace = false;
      tmp[j] = str[i];
      j++;
    }

  }
  tmp[j] = '\0';
  //sprintf(str,"");
  memcpy(str, tmp, READBUFFER);  
  
}



int StartsWith(char strIn[], char strRef[])
{
  int i;
  
  for(i=0;i<strlen(strRef);i++)
    {
      if(strIn[i]!=strRef[i])
	return false;
      //printf("%c %c\n",strIn[i],strRef[i]);
    }
  return true;
}



int FindAtomByName(char name[])
{
  int i;

  for(i=0;i<Run.NTypes;i++)
    if(StartsWith(name,Types[i].Name) == true)
      return i;

  return -1;
}
