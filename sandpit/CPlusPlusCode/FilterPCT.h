#pragma once
//
// FilterPCT.cpp
//
#ifndef __FILTERPCT_H__
#define __FILTERPCT_H__


typedef struct tagPCT_VEC {
	int PCT;
	double dProb;
} PCT_VEC;


int GetNumRows(char *InPath);

int GetNumCols(char *InPath);

bool DirectCopyPCTFile(char *InPath);

bool CreateShapeTable(char *InPath, int NumPCTs);

int getPCTIndex(int pct);

long GetGridOffset(double dX, double dY, EsriBinaryHeader *ebh);

#endif // __FILTERPCT_H__