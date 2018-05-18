//
// gdmTypes.h
//
#ifndef __GDMTYPES_H__
#define __GDMTYPES_H__

#if defined _M_X64

typedef long long GDM_INT;
typedef long long GDM_UINT;

#elif defined _WIN32

typedef int GDM_INT;
typedef unsigned int GDM_UINT;

#endif


#endif // __GDMTYPES_H__