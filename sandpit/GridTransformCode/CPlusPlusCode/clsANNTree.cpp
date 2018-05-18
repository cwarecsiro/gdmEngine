//
// clsANNTree.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "myCallback.h"
#include "Message.h"
#include "ParamsW16.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "clsANNTree.h" 


//
// default ctor
//
ANNTree::ANNTree(char *pDataCSVFile)
{
	nRows = GetNumRows(pDataCSVFile);

	//Message(nRows, "nRows");
	if (nRows < 1 )
	{
		Message("pDataCSVFile has NO row data!", "ERROR");
		return;
	}

	nCols = GetNumCols(pDataCSVFile, 0);
	if (nCols < 1 )
	{
		Message("pDataCSVFile has NO column data!", "ERROR");
		return;
	}
	//Message(nCols, "nCols");	
	//Message(nColsToSkip, "nColsToSkip");

	// allocate and populate the ANN lookup table
	pTableData = annAllocPts( nRows, nCols );
	PopulateTable(pDataCSVFile, pTableData, 0, nRows, nCols);
	//PrintLookupTable("c:\\data\\sga_impact_2015\\aaTestPrintout_01.csv");

	nn_idx = new ANNidx[nRows];
	dist = new ANNdist[nRows];
	query_pt = annAllocPt(nCols);

	theTree = new ANNkd_tree ( pTableData, nRows, nCols ); 
}


//
// ctor based on number of leading columns to skip (zero based)
//
ANNTree::ANNTree(char *pDataCSVFile, int nColsToSkip)
{
	nRows = GetNumRows(pDataCSVFile);
	//Message(nRows, "nRows");
	if (nRows < 1 )
	{
		Message("pDataCSVFile has NO row data!", "ERROR");
		return;
	}

	nCols = GetNumCols(pDataCSVFile, nColsToSkip);
	if (nCols < 1 )
	{
		Message("pDataCSVFile has NO column data!", "ERROR");
		return;
	}
	//Message(nCols, "nCols");	

	// allocate and populate the ANN lookup table
	pTableData = annAllocPts(nRows, nCols);
	PopulateTable(pDataCSVFile, pTableData, nColsToSkip, nRows, nCols);
	//PrintLookupTable("c:\\data\\sga_impact_2015\\aaTestPrintout_01.csv");

	nn_idx = new ANNidx[nRows];
	dist = new ANNdist[nRows];
	query_pt = annAllocPt(nCols);
	theTree = new ANNkd_tree ( pTableData, nRows, nCols ); 

	/*const int arraySize = 5;
	int testIndices[arraySize] = {5,10,20,50,75};
	char filepath[BUFFLEN];
	for (int i=0; i<arraySize; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			query_pt[j] = pTableData[testIndices[i]][j];
		}
		Search();

		if (testIndices[i] < 10)
			sprintf(filepath, "c:\\data\\sga_impact_2015\\aaTestPrintout_0%d.csv", testIndices[i]); 
		else
			sprintf(filepath, "c:\\data\\sga_impact_2015\\aaTestPrintout_%d.csv", testIndices[i]); 

		PrintIndexAndDistVectors(filepath);
	}*/
}


//
// dtor
//
ANNTree::~ANNTree()
{
	if (theTree) delete theTree;
}


//
// perform the ANN search
//
bool 
ANNTree::Search()
{
	GetTree()->annkSearch( GetQueryPt(), GetNumRows(), GetIndices(), GetDist(), 0 );
	return(true);
}


//
// counts rows but assume presence of and skip header row
//
int 
ANNTree::GetNumRows(char *pDataCSVFile)
{
	// open the site file
	FILE *fp = fopen( pDataCSVFile, "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	// determine the number of records in the site file
	int nItems = 0;
	while(1)
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		++nItems;
	}

	// close file
	fclose( fp );
	return(nItems);
}


//
// Get the number of columns in the demand site file 
//
int 
ANNTree::GetNumCols(char *pDataCSVFile, int nColsToSkip)
{
	// open the site file
	FILE *fp = fopen( pDataCSVFile, "r+t" );

	// get the header and close file
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );
	fclose( fp );

	int nColumns = 0;

	// the case where we use all the columns
	if ( nColsToSkip < 1 )
	{
		char *p = strtok(buff, ",\n");
		if ( p ) ++nColumns;
		do 
		{
			p = strtok(NULL, ",\n");
			if (p) 
				++nColumns;

		} while(p);
	}

	else // we need to skip nColsToSkip leading columns
	{
		int nSkipCols = 0;
		char *p = strtok(buff, ",\n");
		if ( p )
		{			
			++nSkipCols;
		}
		
		if ( nColsToSkip > 1 )
		{
			do 
			{
				if (nSkipCols >= nColsToSkip)
					break;

				p = strtok(NULL, ",\n");
				if (p) ++nSkipCols;
			} while(p);
		}

		// now we've skipped the desired number of columns, count what's left...
		do 
		{
			p = strtok(NULL, ",\n");
		
			if (p) 
				++nColumns;

		} while(p);
	}

	return(nColumns);
}


//
// Populate the ANN lookup table fornt e CSV file
//
void 
ANNTree::PopulateTable(char *pDataCSVFile, ANNpointArray pTableData, int nColsToSkip, int nRows, int nCols)
{
	// open the site file
	FILE *fp = fopen( pDataCSVFile, "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	for ( int i=0; i<nRows; i++ )
	{
		fgets( buff, BUFFLEN, fp );

		if (nColsToSkip < 1)
		{
			char *p = strtok(buff, ",\n");
			if ( p ) 
			{
				pTableData[i][0] = atof(p);
			}

			for ( int j=1; j<nCols; j++ )
			{
				p = strtok(NULL, ",\n");
				pTableData[i][j] = atof(p);
			}
		}

		else
		{
			char *p = strtok(buff, ",\n");
			for ( int j=1; j<nColsToSkip; j++ )
			{
				p = strtok(NULL, ",\n");
			}

			for ( int j=0; j<nCols; j++ )
			{
				p = strtok(NULL, ",\n");
				pTableData[i][j] = atof(p);
			}
		}
	}

	fclose( fp );
}


//
// Print the table to csv file
//
void 
ANNTree::PrintLookupTable(char *pOutFile)
{
	FILE *fp = fopen(pOutFile, "w+t");
	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			fprintf(fp, "%lf", pTableData[i][j]);

			if (j < nCols-1 )
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
}


//
// Print the indices and distances to a csv file
//
void 
ANNTree::PrintIndexAndDistVectors(char *pOutFile)
{
	FILE *fp = fopen(pOutFile, "w+t");
	fprintf(fp, "Index,Distance,Dissimilarity\n");
	for ( int i=0; i<nRows; i++ )
	{
		fprintf(fp, "%d,%lf,%lf\n", nn_idx[i], dist[i], 1-exp(-dist[i]));
	}
	fclose(fp);
}
