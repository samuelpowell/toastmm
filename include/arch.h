// =======================================================================
// Architecture-dependent features
// =======================================================================

#ifndef __TOASTARCH_H
#define __TOASTARCH_H

#if defined(_WIN32) || defined(_WIN64)

#include <process.h>	// getpid
#include <direct.h>		// getcwd
#include <float.h>		// rng
#include <stdlib.h>

#define __func__ __FUNCTION__

// MSVC string workarounds
#define strcasecmp(str1,str2) _stricmp((str1),(str2))
#define strncasecmp(str1,str2,n) _strnicmp((str1),(str2),(n))

// Avoid MSVC deprecation warnings
#define getpid() _getpid()
#define getcwd(buffer,maxlen) _getcwd(buffer,maxlen)
#define unlink(fname) _unlink(fname)

// MSVC replacement for non-standard RNG
inline double drand48(void) { return (double)rand()/(double)RAND_MAX; }

#endif // WINxx

#endif // !__TOASTARCH_H
