#ifndef READBUFFER



#define true 1
#define false 0
#define READBUFFER 200



typedef struct
{
  char Name[5];
  float Mass;
  float Charge;
}AtomType;


typedef struct
{
  int ID;
  int TypeID;
  int ClusterID;
  float poscore[3];
  float Force[3];
}Atom;

typedef struct
{
  int NTypes;
  int NAtoms;
  int NCls;

  unsigned short NTipAtoms;

  int DeviceID;

  int MDimages;
  int BperIMG,FperIMG;

  int MD_SaveXYZ, MD_XYZfreq;
  int MD_equiloop;

  float MD_Step;
  int MD_NSteps, MD_ESteps;
  int MD_Stats;

  int MD_reload;
  int MD_reload_num;
  char MD_reloadfile[READBUFFER];
  int MD_reinit;

  int ChainReset;
  int TipDyn;
  int   TipOscAxis;
  float TipDynPars[5];
  float TipDynPar2[5];

  short TipBound; //number of tip atoms
  int TipCls;
  int TipClsIdx[10];
  float TipXY[3];
  float TipXY0[3];

  int ThermoTipCls, ThermoSurfCls;
  int ThermoTip[5], ThermoSurf[5]; //indexes of tip and surface thermostat clusters
  float TempTip, TempSurf;
  float TCoupTip,TCoupSurf;  //berendsen coupling strengths

}RunParams;

typedef struct
{
  char Name[20];
  short hfix;
  float spring;  
}ClusterDef;


typedef struct
{
  unsigned int Type;
  float Value[6];
  
  
}PairPot;

//Atom *AtomSystem;


//GLOBAL VARIABLES
//int *SMap;

unsigned int sSize;
RunParams Run;
AtomType *Types;
Atom *tSystem, *System;
//double *Dists;  //array with distances
ClusterDef *Clusters;
PairPot *PairPots;




//FUNCTIONS//

int main(int argc, char *argv[]);
int Map(int i, int j);

int DeAllocate();

int FindAtomByName(char name[]);


/*
//GPU Functions//
int ProbeDevice(int device);
int AllocateGPU();
int FreeGPU();

int GetForces(int CopyBack);
int DoOptStep(int CopyBack);

*/




//Input Readers
int CheckInput(int argc, char *argv[]);



//Output Functions
int DBG_Atoms();




#endif
