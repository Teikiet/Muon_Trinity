/*----------------------------------------------------------------------------
*
* Copyright (C) 2021 - 2024 Antonio Augusto Alves Junior
*
* This file is part of RandomIterator.

* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* 3. Neither the name of the copyright holder nor the names of its contributors
* may be used to endorse or promote products derived from this software without
* specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
* STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
* OF THE POSSIBILITY OF SUCH DAMAGE.
*
* ----------------------------------------------------------------------------*/

/*
 * RandomIterator.h
 *
 *  Created on: 21/02/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

//  This is the only RandomIterator header that is guaranteed to
//  change with every RandomIterator release.
//
//  RandomIterator_VERSION % 100 is the sub-minor version
//  RandomIterator_VERSION / 100 % 1000 is the minor version
//  RandomIterator_VERSION / 100000 is the major version
//
//  Because this header does not #include <RandomIterator/detail/Config.h>,
//  it is the only RandomIterator header that does not cause
//  RandomIterator_HOST_SYSTEM and RandomIterator_DEVICE_SYSTEM to be defined.

/*! \def RandomIterator_VERSION
 *  \brief The preprocessor macro \p RandomIterator_VERSION encodes the version
 *         number of the RandomIterator.
 *
 *         <tt>RandomIterator_VERSION % 100</tt> is the patch version.
 *         <tt>RandomIterator_VERSION / 100 % 1000</tt> is the feature version.
 *         <tt>RandomIterator_VERSION / 100000</tt> is the major version.
 */
#define RandomIterator_VERSION 100001


/*! \def RandomIterator_MAJOR_VERSION
 *  \brief The preprocessor macro \p RandomIterator_MAJOR_VERSION encodes the
 *         major version number of RandomIterator.
 */
#define RandomIterator_MAJOR_VERSION     (RandomIterator_VERSION / 100000)

/*! \def RandomIterator_MINOR_VERSION
 *  \brief The preprocessor macro \p RandomIterator_MINOR_VERSION encodes the
 *         minor version number of RandomIterator.
 */
#define RandomIterator_MINOR_VERSION     (RandomIterator_VERSION / 100 % 1000)

/*! \def RandomIterator_PATCH_NUMBER
 *  \brief The preprocessor macro \p RandomIterator_PATCH_NUMBER encodes the
 *         patch number of the RandomIterator library.
 */
#define RandomIterator_PATCH_NUMBER  0



// Declare these namespaces here for the purpose of Doxygenating them

/*!
 *  \brief \p random_iterator is the top-level namespace which contains all RandomIterator
 *         functions and types.
 */
namespace random_iterator{ }


