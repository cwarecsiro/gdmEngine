//
// ConvertToSxS.h
//
#ifndef __CONVERTTOSXS_H__
#define __CONVERTTOSXS_H__

#include "myCallback.h"

extern "C" 
{ 

	//
	// Convert a Comma delimited textfile in the form X,Y,Species_Code
	// to a Site By Species File in the form X,Y,Sp1,Sp2,...,SpN where
	// Sp1..SpN are presence/absence values (1/0)
	//
	_declspec(dllexport) bool SitePlusSpeciesToSxS(char *pInput, char *pOutput, FPTR fptr);


	//
	// Convert a Comma delimited textfile in the form X,Y,Species_Code
	// to a Site By Species File in the form X,Y,Sp1,Sp2,...,SpN where
	// Sp1..SpN are abundance values (0..255)
	//
	_declspec(dllexport) bool SitePlusSpeciesAbundanceToSxS(char *pInput, char *pOutput, FPTR fptr);


	//
	// Get the Minimum and Maximum counts from a Site By Species Table
	//
	_declspec(dllexport) bool GetSxSMinMax(char *pInput, int *pMin, int *pMax, FPTR fptr);

}

#endif // __CONVERTTOSXS_H__