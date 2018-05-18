//
// TestSubSample.h
//
#ifndef __TESTSUBSAMPLE_H__
#define __TESTSUBSAMPLE_H__

extern "C" 
{
	//
	// Test the subsampling code
	//
	_declspec(dllexport) bool TestSubsample_01(char *sPath, int Count, int Min, int Max);
}



#endif // __TESTSUBSAMPLE_H__