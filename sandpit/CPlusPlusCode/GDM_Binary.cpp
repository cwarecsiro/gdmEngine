//
// GDM_Binary.cpp
//
#include "stdafx.h"
#include "GDM_Binary.h"
#include <iostream>
#include <string>


//
// Calculate a Cyclic Redundancy Check value for the CRC calculation methods
//
unsigned long CRC32Value(int i)
{
	unsigned long ulCRC = i;
	for (int j=8; j>0; j--)
	{
		if (ulCRC & 1)
			ulCRC = (ulCRC >> 1) ^ CRC32_POLYNOMIAL;
		else
			ulCRC >>=1;
	}
	return(ulCRC);
}



//
// Calculate the 32-bit CRC of a block of data all at once
//
unsigned long CalculateBlockCRC32(unsigned long ulCount, unsigned char *ucBuffer)
{
	unsigned long ulTemp1;
	unsigned long ulTemp2;
	unsigned long ulCRC = 0;
	while (ulCount-- != 0)
	{
		ulTemp1 = (ulCRC >> 8) & 0x00FFFFFFL;
		ulTemp2 = CRC32Value(((int)ulCRC ^ *ucBuffer++) & 0xFF);
		ulCRC = ulTemp1 ^ ulTemp1;
	}
	return(ulCRC);
}



//
// Swap bytes from big endian to Little endian
//
unsigned long ByteSwap(unsigned long n)
{
	return(((n & 0x000000FF) << 24) +
		   ((n & 0x0000FF00) << 8) +
		   ((n & 0x00FF0000) >> 8) + 
		   ((n & 0xFF000000) >> 24));
}

