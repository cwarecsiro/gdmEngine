//
//
// clsEsriBinaryHeader.cpp
//
//
#include "stdafx.h"

#include "clsEsriBinaryHeader.h"
#include "Message.h"

//
// ctor - Opens the .hdr file and extract the respective grid parameters
//        If the lower left cell locator is in the gridcell center, then adjust to lower left cell corner
//        Set the byte order enum to 0 or 1 depending on the byte order of the binary data in the .flt file
//
EsriBinaryHeader::EsriBinaryHeader(char *hdrFile)
{
	if ( false == Init(hdrFile) )
	{
		Message("<EsriBinaryHeader ctor> Cannot do Init()", "ERROR" );
		return;
	}

	if ( false == SetValuesFromFile() )
	{
		Message("<EsriBinaryHeader ctor> Cannot SetValuesFromFile()", "ERROR" );
		return;
	}
}; 


//
// dtor
//
EsriBinaryHeader::~EsriBinaryHeader()
{		
	if ( FilePath) delete[] FilePath;
};



//
// Initialse the private internals
//
bool EsriBinaryHeader::Init(char *hdrFile)
{
	FilePath = new char [ESRI_MAX_BUFF_SIZE];
	if ( NULL == FilePath )
	{
		Message("<EsriBinaryHeader::Init> Cannot allocate FilePath", "ERROR" );
		return(false);
	}
	strcpy(FilePath, hdrFile);

	NumRows = 0;
	NumCols = 0;
	CellSize = 0.0;
	XllCorner = 0.0;
	YllCorner = 0.0;
	NoDataValue = 0.0F;
	ByteOrder = -1;
	return(true);
}


//
// Display the private internal grid values
//
void EsriBinaryHeader::Display()
{
	char *p = new char [ESRI_MAX_BUFF_SIZE];
	if ( NULL == p )
	{
		Message("<EsriBinaryHeader::Display> Cannot allocate p", "ERROR" );
		return;
	}

	sprintf( p, 
		     "FilePath: %s\n\nNumCols: %d\nNumRows: %d\nXllCorner: %lf\nYllCorner: %lf\nCellSize: %lf\nNoData: %f\nByteOrder: %s",
			 FilePath,GetNumCols(), GetNumRows(), GetXllCorner(), GetYllCorner(), GetCellSize(), GetNoDataValue(), 
			 (ByteOrder==0)?"LSBFIRST":(ByteOrder==0)?"MSBFIRST":"UNINITIALISED");

	Message(p, "Esri Binary Grid File Header" );
	if ( p ) delete[] p;
}


//
// The following methods extract the respective fields and set values from the header file
//
bool EsriBinaryHeader::SetValuesFromFile()
{
	if ( false == SetNumRowsFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract nrows from header file", "ERROR" );
		return(false);
	}

	if ( false == SetNumColsFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract ncols from header file", "ERROR" );
		return(false);
	}

	if ( false == SetCellSizeFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract cellsize from header file", "ERROR" );
		return(false);
	}

	if ( false == SetXllCornerFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract xllcorner from header file", "ERROR" );
		return(false);
	}

	if ( false == SetYllCornerFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract yllcorner from header file", "ERROR" );
		return(false);
	}	

	if ( false == SetNoDataFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract nodata from header file", "ERROR" );
		return(false);
	}

	if ( false == SetByteOrderFromFile() )
	{
		Message("<EsriBinaryHeader::SetValuesFromFile> Cannot extract byteorder from header file", "ERROR" );
		return(false);
	}

	return(true);
}


bool EsriBinaryHeader::SetNumRowsFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetNumRowsFromFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive nrows
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetNumRowsFromFile()> At EOF without locating nrows", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "nrows") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	NumRows = atoi(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetNumColsFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetNumColsFromFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive ncols
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetNumColsFromFile()> At EOF without locating nCols", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "ncols") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	NumCols = atoi(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetXllCornerFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetXllCornerFomFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive xcorner
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);

			// try accessing the xllcenter field and making adjustment for the xllcorner value
			if ( false == SetXllCornerFromXllCenter() )
			{
				Message("<EsriBinaryHeader::SetXllCornerFomFile()> At EOF without locating XllCorner ot XllCenter", "ERROR" );			
				return(false);
			}
			else
			{
				// we have successfully acessed the xllcenter and modified it for the xllcorner value
				return(true);
			}
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "xllcorner") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	XllCorner = atof(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetXllCornerFromXllCenter()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetXllCornerFromXllCenter()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive xcorner
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetXllCornerFromXllCenter()> At EOF without locating XllCenter", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "xllcenter") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	XllCorner = atof(p) - (GetCellSize() / 2.0);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetYllCornerFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetYllCornerFomFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive yllcorner
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);

			// try accessing the xllcenter field and making adjustment for the xllcorner value
			if ( false == SetYllCornerFromYllCenter() )
			{
				Message("<EsriBinaryHeader::SetYllCornerFomFile()> At EOF without locating YllCorner ot YllCenter", "ERROR" );			
				return(false);
			}
			else
			{
				// we have successfully acessed the yllcenter and modified it for the yllcorner value
				return(true);
			}
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "yllcorner") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	YllCorner = atof(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetYllCornerFromYllCenter()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetYllCornerFromYllCenter()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive xcorner
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetYllCornerFromYllCenter()> At EOF without locating YllCenter", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "yllcenter") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	YllCorner = atof(p) - (GetCellSize() / 2.0);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetCellSizeFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetCellSizeFomFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive cellsize
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetCellSizeFomFile()> At EOF without locating CellSize", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "cellsize") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	CellSize = (float)atof(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetNoDataFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetNoDataFromFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive nodata
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetNoDataFromFile()> At EOF without locating NoData", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "nodata_value") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	NoDataValue = (float)atof(p);

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


bool EsriBinaryHeader::SetByteOrderFromFile()
{
	char *p;
	char *buff = new char[ESRI_MAX_BUFF_SIZE];

	FILE *fp = fopen(FilePath,"r+t");
	if ( NULL == fp )
	{
		Message("<EsriBinaryHeader::SetByteOrderFromFile()> Cannot open FilePath", "ERROR" );
		if ( buff) delete[] buff;
		return(false);
	}

	// locate the line containing case insensitive byteorder
	while(1)
	{
		if ( NULL == fgets(buff, ESRI_MAX_BUFF_SIZE, fp) )
		{
			Message("<EsriBinaryHeader::SetByteOrderFromFile()> At EOF without locating ByteOrder", "ERROR" );
			if ( buff) delete[] buff;
			if ( fp ) fclose(fp);
			return(false);
		}

		p = strtok(buff," \t\n");
		if ( 0 == _stricmp(p, "byteorder") )
		{
			break;
		}
	}
	
	// extract the value
	do 
	{
		p = strtok( NULL, " \t\n");
	} while (0 == strcmp(p,""));
	if ( 0 == _stricmp(p, "LSBFIRST") )
		ByteOrder = 0;
	else if ( 0 == _stricmp(p, "MSBFIRST") )
		ByteOrder = 1;
	else
	{
		Message("<EsriBinaryHeader::SetByteOrderFromFile> Cannot decode byteorder", "ERROR" );
		if ( fp ) fclose(fp);
		if ( buff) delete[] buff;
		return(false);
	}

	// cleanup
	if ( fp ) fclose(fp);
	if ( buff) delete[] buff;
	return(true);
}


//
// The following public access methods return the binary grid metrics to the class owner
//
int EsriBinaryHeader::GetNumRows()
{
	return(NumRows);	
};


int EsriBinaryHeader::GetNumCols()
{
	return(NumCols);		
};


double EsriBinaryHeader::GetCellSize()
{
	return(CellSize);		
};


double EsriBinaryHeader::GetXllCorner()
{
	return(XllCorner);		
};


double EsriBinaryHeader::GetYllCorner()
{
	return(YllCorner);		
};



float EsriBinaryHeader::GetNoDataValue()
{
	return(NoDataValue);		
};


//bool EsriBinaryHeader::CopyTo( char *OutPath )
//{
//	FILE *fp = fopen( OutPath, "w+t" );
//	if ( NULL == fp )
//	{
//		Message("<EsriBinaryHeader::CopyTo> - Cannot open OutPath", "ERROR");
//		return(false);
//	}
//	fprintf( fp, "ncols         %d\n", this->NumCols );
//	fprintf( fp, "nrows         %d\n", this->NumRows );
//	fprintf( fp, "xllcorner     %lf\n", this->XllCorner );
//	fprintf( fp, "yllcorner     %lf\n", this->YllCorner );
//	fprintf( fp, "cellsize      %lf\n", this->CellSize );
//	fprintf( fp, "NODATA_value  %1.0f\n", this->NoDataValue );
//	if ( this->ByteOrder == LSBFIRST )
//	{
//		fprintf( fp, "byteorder     LSBFIRST\n" );
//	}
//	else if ( this->ByteOrder == MSBFIRST )
//	{
//		fprintf( fp, "byteorder     MSBFIRST\n" );
//	}
//	else
//	{
//		Message("<EsriBinaryHeader::CopyTo()> - Invalid Esri Byte Order", "ERROR");
//		fclose(fp);
//		return(false);
//	}
//	fclose(fp);
//	return(true);
//}

bool EsriBinaryHeader::CopyTo(char *OutPath)
{
	FILE *fpOut = fopen(OutPath, "w+t");
	if (NULL == fpOut)
	{
		Message("<EsriBinaryHeader::CopyTo> - Cannot open OutPath", "ERROR");
		return(false);
	}

	FILE *fpIn = fopen(FilePath, "r+t");
	if (NULL == fpIn)
	{
		Message("<EsriBinaryHeader::CopyTo> - Cannot open FilePath for Input", "ERROR");
		return(false);
	}

	char buff[1024];
	while (1)
	{
		if (NULL == fgets(buff, 1024, fpIn))
		{
			break;
		}

		fprintf(fpOut, "%s", buff);
	}

	fclose(fpIn);
	fclose(fpOut);
	return(true);
}


//
// ROI properties
//
double EsriBinaryHeader::GetMinX()
{
	return(XllCorner);		
}


double EsriBinaryHeader::GetMaxX()
{
	return((NumCols * CellSize) + XllCorner);
}


double EsriBinaryHeader::GetMinY()
{
	return(YllCorner);		
}


double EsriBinaryHeader::GetMaxY()
{
	return((NumRows * CellSize) + YllCorner);
}
