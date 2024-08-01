#define CAT 1
#define ORD 2

#define BIRTH 1
#define DEATH 2
#define SWAP 3
#define CHANGE 4

#define NORMAL_REG		  1
#define LOGIT_REG		  2
#define POISSON_REG		  3

extern double pBD;
extern double pSwap;
extern double pChange;
extern double aa;
extern double bb;
extern int NumX;
extern int *VarType;
extern int NumY;
extern int NumDe;
extern int NumWi;
extern int NumV;

extern int NumObs;
extern double **XDat;
extern double **YDat;
extern double **DeDat;
extern double **WiDat;
extern double **VDat;
extern double* YDat1;
extern double* DeDat1;
extern double* WiDat1;
extern double* VDat1;
extern double **XDatR;
extern int NumXR;
extern double* weights;

extern int *RuleNum;
extern double **RuleMat;

extern int *Ivec; //global vector defined to be 1,2,3,...NumObs in main

#include "CPriParams.h"
#include "EndNodeModel.h"

extern CPriParams PriParams;
extern EndNodeModel* endNodeModel;



