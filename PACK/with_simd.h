// control usage of Intel SIMD intrinsics
// by default SIMD intrinsics and assembler code are used on x86_64 and aarch64
// use -DNO_SIMD to refrain from using SIMD intrinsics on x86_64 architectures
#if ! defined(NO_SIMD)
#define WITH_SIMD
// #pragma message("NOTE: using SIMD intrinsics, use -DNO_SIMD to use pure C code")
#warning "NOTE: using Intel SIMD intrinsics, use -DNO_SIMD to use pure C code"
#endif
