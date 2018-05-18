//
//
// clsEsriBinaryHeader.h
//
//

#ifndef	__CLSESRIBINARYHEADER_H__
#define	__CLSESRIBINARYHEADER_H__

#define ESRI_MAX_BUFF_SIZE 500
#define LSBFIRST 0
#define MSBFIRST 1

class EsriBinaryHeader
{
public:

	EsriBinaryHeader(char *hdrFile); 

	~EsriBinaryHeader();

	int GetNumRows();

	int GetNumCols();

	double GetCellSize();

	double GetXllCorner();

	double GetYllCorner();

	float GetNoDataValue();

	int GetByteOrder();

	void Display();

	bool CopyTo( char *OutPath );

	double GetMinX();

	double GetMaxX();

	double GetMinY();

	double GetMaxY();


private:

	bool Init(char *hdrFile);

	bool SetValuesFromFile();

	char *FilePath;

	int NumRows;

	int NumCols;

	double CellSize;

	double XllCorner;

	double YllCorner;

	float NoDataValue;

	int ByteOrder;

	bool SetNumRowsFromFile();
	bool SetNumColsFromFile();
	bool SetXllCornerFromFile();
	bool SetYllCornerFromFile();
	bool SetXllCornerFromXllCenter();
	bool SetYllCornerFromYllCenter();
	bool SetCellSizeFromFile();
	bool SetNoDataFromFile();
	bool SetByteOrderFromFile();
};


#endif //__CLSESRIBINARYHEADER_H__