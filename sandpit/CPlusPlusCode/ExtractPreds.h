//
// ExtractPreds.h
//
#ifndef __EXTRACTPREDS_H__
#define __EXTRACTPREDS_H__

double **ExtractLookupTable(char *sPath, bool DoBatch, int *pSize);

bool ValidateLookupFile(char *pPath, bool DoBatch, int *pSize);


#endif // __EXTRACTPREDS_H__