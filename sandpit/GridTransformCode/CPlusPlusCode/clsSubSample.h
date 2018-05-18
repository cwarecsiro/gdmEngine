//
// clsSubSample.h
//
#ifndef __CLSSUBSAMPLE_H__
#define __CLSSUBSAMPLE_H__

class SubSample
{
	 unsigned *pNumbers;
	 int SampleSize;
	 int RangeMin;
	 int RangeMax;
	 bool Valid;

public: 
	SubSample(int nSampleSize, int nRangeMin, int nRangeMax);

	bool IsValid() {return(Valid);}
	int GetSampleSize() {return(SampleSize);}
	int GetRangeMin() {return(RangeMin);}
	int GetRangeMax() {return(RangeMax);}
	unsigned GetNumberAt(int Index) {return(pNumbers[Index]);}
	void Sort();
	
	void Stats();
	void Print();
	void WriteToText(char *sPath);
	~SubSample();

};



#endif // __CLSSUBSAMPLE_H__