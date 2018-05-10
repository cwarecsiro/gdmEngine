//
// clsMemorySxSDouble.h
//
#ifndef __CLSMEMORYSXSDOUBLE_H__	
#define __CLSMEMORYSXSDOUBLE_H__


class MemorySxSDouble
{
	double *pX;
	double *pY;
	double **ppSxS;
	int numRows;
	int numDataCols;
	bool Valid;

	public: 
		MemorySxSDouble(char *pPath);

		int GetRowCount() {return(numRows);}
		int GetColCount() {return(numDataCols);}
		double **GetSxS() {return(ppSxS);}
		double *GetX() {return(pX);}
		double *GetY() {return(pY);}
		bool IsValid() {return(Valid);}

		~MemorySxSDouble();

	private:

		int CountDataColumns(char *pPath);
		int CountRows(char *pPath);
		bool PopulateFromSxS(char *pPath);
};



#endif // __CLSMEMORYSXSDOUBLE_H__	