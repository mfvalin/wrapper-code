// control usage of Intel SIMD intrinsics
// by default SIMD intrinsics and assembler code are used on x86_64 and aarch64
// use -DNO_SIMD to refrain from using SIMD intrinsics on x86_64 architectures
#if ! defined(WITH_SIMD_CONTROL_DONE)
#define WITH_SIMD_CONTROL_DONE

#if ! defined(NO_SIMD)
#define WITH_SIMD

#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#endif
#if defined(__x86_64__) && defined(__SSE2__)
#include <emmintrin.h>
#endif

#if ! defined(QUIET_SIMD)
// #pragma message("NOTE: using SIMD intrinsics, use -DNO_SIMD to use pure C code")
#warning "NOTE: using Intel SIMD intrinsics, use -DNO_SIMD to use pure C code"
#endif

#endif

#endif
