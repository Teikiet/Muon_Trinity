/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef __sunprofeatures_dot_hpp
#define __sunprofeatures_dot_hpp

#ifndef RANDOM_ITERATOR_R123_STATIC_INLINE
#define RANDOM_ITERATOR_R123_STATIC_INLINE static inline
#endif

#ifndef RANDOM_ITERATOR_R123_FORCE_INLINE
#define RANDOM_ITERATOR_R123_FORCE_INLINE(decl) decl
#endif

#ifndef RANDOM_ITERATOR_R123_CUDA_DEVICE
#define RANDOM_ITERATOR_R123_CUDA_DEVICE
#endif

#ifndef RANDOM_ITERATOR_R123_ASSERT
#include <assert.h>
#define RANDOM_ITERATOR_R123_ASSERT(x) assert(x)
#endif

#ifndef RANDOM_ITERATOR_R123_BUILTIN_EXPECT
#define RANDOM_ITERATOR_R123_BUILTIN_EXPECT(expr,likely) expr
#endif

// The basic idiom is:
// #ifndef RANDOM_ITERATOR_R123_SOMETHING
// #if some condition
// #define RANDOM_ITERATOR_R123_SOMETHING 1
// #else
// #define RANDOM_ITERATOR_R123_SOMETHING 0
// #endif
// #endif
// This idiom allows an external user to override any decision
// in this file with a command-line -DRANDOM_ITERATOR_R123_SOMETHING=1 or -DRANDOM_ITERATOR_R123_SOMETHINE=0

// An alternative idiom is:
// #ifndef RANDOM_ITERATOR_R123_SOMETHING
// #define RANDOM_ITERATOR_R123_SOMETHING (some boolean expression)
// #endif
// where the boolean expression might contain previously-defined RANDOM_ITERATOR_R123_SOMETHING_ELSE
// pp-symbols.

#ifndef RANDOM_ITERATOR_R123_USE_AES_NI
#define RANDOM_ITERATOR_R123_USE_AES_NI 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_SSE4_2
#define RANDOM_ITERATOR_R123_USE_SSE4_2 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_SSE4_1
#define RANDOM_ITERATOR_R123_USE_SSE4_1 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_SSE
#define RANDOM_ITERATOR_R123_USE_SSE 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_AES_OPENSSL
#define RANDOM_ITERATOR_R123_USE_AES_OPENSSL 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_GNU_UINT128
#define RANDOM_ITERATOR_R123_USE_GNU_UINT128 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_ASM_GNU
#define RANDOM_ITERATOR_R123_USE_ASM_GNU 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_CPUID_MSVC
#define RANDOM_ITERATOR_R123_USE_CPUID_MSVC 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_X86INTRIN_H
#define RANDOM_ITERATOR_R123_USE_X86INTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_IA32INTRIN_H
#define RANDOM_ITERATOR_R123_USE_IA32INTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_XMMINTRIN_H
#define RANDOM_ITERATOR_R123_USE_XMMINTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_EMMINTRIN_H
#define RANDOM_ITERATOR_R123_USE_EMMINTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_SMMINTRIN_H
#define RANDOM_ITERATOR_R123_USE_SMMINTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_WMMINTRIN_H
#define RANDOM_ITERATOR_R123_USE_WMMINTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_INTRIN_H
#define RANDOM_ITERATOR_R123_USE_INTRIN_H 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO16_ASM
#define RANDOM_ITERATOR_R123_USE_MULHILO16_ASM 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO32_ASM
#define RANDOM_ITERATOR_R123_USE_MULHILO32_ASM 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO64_ASM
#define RANDOM_ITERATOR_R123_USE_MULHILO64_ASM 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO64_MSVC_INTRIN
#define RANDOM_ITERATOR_R123_USE_MULHILO64_MSVC_INTRIN 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO64_CUDA_INTRIN
#define RANDOM_ITERATOR_R123_USE_MULHILO64_CUDA_INTRIN 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_MULHILO64_OPENCL_INTRIN
#define RANDOM_ITERATOR_R123_USE_MULHILO64_OPENCL_INTRIN 0
#endif

#ifndef RANDOM_ITERATOR_R123_USE_PHILOX_64BIT
#define RANDOM_ITERATOR_R123_USE_PHILOX_64BIT 0
#endif

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>
#ifndef UINT64_C
#error UINT64_C not defined.  You must define __STDC_CONSTANT_MACROS before you #include <stdint.h>
#endif

// If you add something, it must go in all the other XXfeatures.hpp
// and in ../ut_features.cpp
#endif
