//
// clsMemorySxS.h
//
#ifndef __CLSMEMORYSXS_H__	
#define __CLSMEMORYSXS_H__

class MemorySxS
{
	double *pX;
	double *pY;
	int **ppSxS;
	int numRows;
	int numDataCols;
	bool Valid;

	public: 
		MemorySxS(char *pPath);

		int GetRowCount() {return(numRows);}
		int GetColCount() {return(numDataCols);}
		int **GetSxS() {return(ppSxS);}
		double *GetX() {return(pX);}
		double *GetY() {return(pY);}
		bool IsValid() {return(Valid);}

		~MemorySxS();

	private:

		int CountDataColumns(char *pPath);
		int CountRows(char *pPath);
		bool PopulateFromSxS(char *pPath);
};


#endif // __CLSMEMORYSXS_H__