//
// gdm2ibra.h
//
#ifndef __GDM2IBRA_H__
#define __GDM2IBRA_H__

#include <stdio.h>
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"
#include "myCallback.h"

int GetIbraLookupMaxIndex(char *IbraPath);

char **GetLookupStrings(char *IbraPath, int nMaxIbra);

void WriteOutput(char *OutPath, 
	             int nCols, 
				 int **ppData, 
				 char **ppLookup, 
				 int *pTallies, 
				 bool CountCells, 
				 double CellSize, 
				 FPTR fptr);


#endif // __GDM2IBRA_H__


