//
// MATLAB Compiler: 8.0 (R2020a)
// Date: Mon Jul 26 18:51:09 2021
// Arguments:
// "-B""macro_default""-W""cpplib:libMatlabFE,all""-T""link:lib""-d""/home/qs/Do
// kumente/Sim3_0_testFunc/libMatlabFE/for_testing""-v""/home/qs/Dokumente/Sim3_
// 0_testFunc/particleThermalComputation.m"
//

#ifndef libMatlabFE_h
#define libMatlabFE_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libMatlabFE_C_API 
#define LIB_libMatlabFE_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_libMatlabFE_C_API 
bool MW_CALL_CONV libMatlabFEInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libMatlabFE_C_API 
bool MW_CALL_CONV libMatlabFEInitialize(void);

extern LIB_libMatlabFE_C_API 
void MW_CALL_CONV libMatlabFETerminate(void);

extern LIB_libMatlabFE_C_API 
void MW_CALL_CONV libMatlabFEPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_libMatlabFE_C_API 
bool MW_CALL_CONV mlxParticleThermalComputation(int nlhs, mxArray *plhs[], int nrhs, 
                                                mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_libMatlabFE
#define PUBLIC_libMatlabFE_CPP_API __declspec(dllexport)
#else
#define PUBLIC_libMatlabFE_CPP_API __declspec(dllimport)
#endif

#define LIB_libMatlabFE_CPP_API PUBLIC_libMatlabFE_CPP_API

#else

#if !defined(LIB_libMatlabFE_CPP_API)
#if defined(LIB_libMatlabFE_C_API)
#define LIB_libMatlabFE_CPP_API LIB_libMatlabFE_C_API
#else
#define LIB_libMatlabFE_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_libMatlabFE_CPP_API void MW_CALL_CONV particleThermalComputation(int nargout, mwArray& TSurface, mwArray& T, const mwArray& nStep, const mwArray& OutStep, const mwArray& dt, const mwArray& HeatFlux, const mwArray& TlastStep);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
