//
// clsSubSample.cpp
//
#include "stdafx.h"
#include "clsSubSample.h"
#include "clsHashtab.h"
#include "Message.h"
#include <time.h>
#include <limits.h>
#include "mtRand.h"


////////////////////////////////////////////////////////////////////////////////////
// String comparison routine for Sort method
//
int uiCompare( const void *arg1, const void *arg2 )
{
	unsigned *p1 = (unsigned *)arg1;
	unsigned *p2 = (unsigned *)arg2;
	if ( *p1 < *p2 )
		return(-1);
	else if ( *p1 > *p2 )
		return(1);
	else
		return(0);
}
////////////////////////////////////////////////////////////////////////////////////


//
// Create a set of randomly generated UNIQUE integers on the half closed interval
// [nRangeMin..nRangeMax) of nSize elements
//
SubSample::SubSample(int nSampleSize, int nRangeMin, int nRangeMax)
{
	SampleSize = nSampleSize;
	RangeMin = nRangeMin;
	RangeMax = nRangeMax;
	

	// make sure that there can be at least nRangeMax integers in the SampleSize
	if ( nSampleSize > abs(nRangeMax-nRangeMin) )
	{
		Message("Requesting MORE random integers that there are UNIQUE numbers", "ERROR");
		Valid = false;
		return;
	}


	//
	// collect the unique random numbers into a Hashtable
	//
	MTRand_int32 *myrand = new MTRand_int32((unsigned)time( NULL ));
	int nItems = 0;	
	Hashtable *hTab = new Hashtable(nSampleSize);
	while(nItems < nSampleSize)
    {
		unsigned int RandVal = myrand->operator()() / 2;
		unsigned int uVal = unsigned((double)RandVal / (INT_MAX) * (nRangeMax - nRangeMin) + nRangeMin);
		if( uVal >= (unsigned)nRangeMax )
			continue;
		char tmp[SIZE_KEY];
		sprintf(tmp, "%d", (int)uVal);
		
		if (!hTab->contains(tmp)) 
		{
			hTab->put(new NODE (tmp,tmp));
			++nItems;
		}
	}


	//
	// assign to the number vector
	//
	pNumbers = new unsigned [nSampleSize];
	if (NULL == pNumbers)
	{
		Message("Cannot allocate pNumbers", "SubSample::SubSample");
		if (hTab) delete hTab;
		Valid = false;
		return;
	}
	int nIndex = 0;
	hTab->initIterator(); 
	while(hTab->hasNext()) 
	{ 
		char qqq[SIZE_KEY];
		hTab->getNextKey(qqq); 
		pNumbers[nIndex++] = atoi(qqq);
	} 
	if (hTab) delete hTab;
	Valid = true;
}



////////////////////////////////////////////////////////////////////////////////////
// Sort the table in ascending order
//
void SubSample::Sort()
{
	qsort( (void *)pNumbers, (size_t)GetSampleSize(), sizeof(unsigned), uiCompare );
}


void SubSample::Stats()
{
	char buff[64];
	sprintf(buff, "Size: %d   Min: %d   Max: %d", GetSampleSize(), GetRangeMin(), GetRangeMax());
	Message(buff, "Stats");
}


void SubSample::Print()
{
	for ( int i=0; i<GetSampleSize(); i++ )
	{
		Message((int)GetNumberAt(i));
	}
}



void SubSample::WriteToText(char *sPath)
{
	FILE *fp = fopen(sPath, "w+t");
	if (NULL == fp)
	{
		Message("Cannot open sPath in SubSample::WriteToText", "ERROR");
		return;
	}

	fprintf( fp, "Index,Value\n" );
	for( int i=0; i<GetSampleSize(); i++ )
	{
		fprintf(fp, "%d,%u\n", i, GetNumberAt(i));
	}
	fclose(fp);
}


SubSample::~SubSample()
{
	if (pNumbers) delete[] pNumbers;
}
