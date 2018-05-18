//
// GDM_Binary.h
//
#ifndef __GDM_BINARY_H__
#define __GDM_BINARY_H__

#define CRC32_POLYNOMIAL 0xEDB88320


//
// Calculate a Cyclic Redundancy Check value for the CRC calculation methods
//
unsigned long CRC32Value(int i);


//
// Calculate the 32-bit CRC of a block of data all at once
//
unsigned long CalculateBlockCRC32(unsigned long ulCount, unsigned char *ucBuffer);


//
// Swap bytes from big endian to Little endian
//
unsigned long ByteSwap(unsigned long n);


#endif // __GDM_BINARY_H__