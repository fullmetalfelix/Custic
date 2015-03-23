// -*- c++ -*-

//#ifndef BLOCKSIZE

#define BLOCKSIZE 64
#define BLOCKFAST 128

//#define MAXSHELLPAR 25
#define READBUFFER 200
#define CHECK 0
#define true 1
#define false 0
#define MAXTYPES 5
#define MAX_PPOT 25 //optimized for 5 species!

#define PI 3.1415926f

#include <cuda.h>
#include <cuda_runtime_api.h>


cudaDeviceProp MyProps;

//standard array sizes
size_t f1_NObj, f1_NAtoms;
size_t f3_NObj, f3_NAtoms;
size_t f4_NObj, f4_NAtoms;
size_t s1_NAtoms;
size_t s3_NAtoms;
size_t i3_NAtoms;
size_t f3_ShellPars;

//System descriptors
float4 *Charges_h, *Charges_d;  // charges positions
float4 *ChargeO_h;              // backup coordinates
float4 *ChargeT_h, *ChargeT_d;  // temporary coordinates
float4 *ChargeINIT;             // atomic positions after the initial setup



short *AtomType_h, *AtomType_d; // integer index for atom ID
short *TipAtoms_h, *TipAtoms_d; //array: element k is 0 if atom k is not part of the tip

float *Mass_h, *Mass_d, *iMass_h, *iMass_d;         // objects masses and inverse masses

float4 *Holders_h, *Holders_d; //positions of holders for fixed atoms (xyz) and spring (w)
float4 *HolderI_d;  //reference position of the holders for the beginning of the loop
float4 *HolderR_d;  //positions of the holders relative to holder[0]

float4 *Speeds_h, *Speeds_d;   // particles speeds (xyz) and mass (w)


float3 *forces_h, *forces_d;   // forces on the objects
//float3 *sforce_h;
//float3 *forceO_h, *ShellForce_d;
//float *Potentials_h, *Potentials_d;


// Temperature measurement
float *Temperature_CPU;
float *Temperature_GPU;       //array with temperature of each image

__constant__ float Tset_d[2];   //temperature setpoints: [0] tip, [1] surface
__constant__ float Tcup_d[2];   //berendsen coupling strengths

__constant__ short TipBound_d;  //number of tip atoms
__constant__ float TipOffset_d; 
// Tip forces stuff
float3 *TipForce_CPU, *TipForce_GPU;        //tip force for each image



//interaction parameters
float3 *ShellParameters_h;      // table with the shell parameters
__constant__ float3 ShellParameters_d[MAX_PPOT];


int FastGrid, FastBlock;
int SlowGrid, SlowBlock;
dim3 TGrid;

__constant__ unsigned short NTipAtoms_d;
__constant__ int NImages_d;
__constant__ int BperIMG_d;
__constant__ int FperIMG_d;

__constant__ int NAtoms_d;
__constant__ int NTypes_d;

short4 ThermoBounds_h;
__constant__ short4 ThermoBounds_d; //index of the tip and surface thermostat regions

__constant__ float MDstep_d,MDQstep_d;
__constant__ float3 TipHrm_d; // paprameters for harmonic tip motion


//__constant__ int LinSize_d;
//__constant__ float alpha;


//__constant__ float2 TypeData[MAXTYPES]; //atom type data (device) [qc qs]


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
  int BperIMG, FperIMG;
 
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
  int TipCls;        //number of tip clusters
  int TipClsIdx[10]; //indexes of the tip clusters
  float TipXY[3];  //where the tip-stickpoint should be after repositioning
  float TipXY0[3]; //where the tip-stickpoint is when the input file is read

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


RunParams Run;
Atom *System;
AtomType *Types;
ClusterDef *Clusters;

//MD calculation variables
int NDegFree;


float HarmonicT0;


//float TipZ, AvgTipZ;
float *TipPosZ_d; //tip position (summed every step for later avg)
float3 TipHldDist; //distance between tip sticker and FirstHolder
float3 TipHld, AvgTipHld;
float3 TipPos, AvgTipPos;





//#endif
