/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_VERSION_H
#define ADIOS_VERSION_H

/* ADIOS Software release version */
#define ADIOS_VERSION "1.9.0"
#define ADIOS_VERSION_MAJOR 1
#define ADIOS_VERSION_MINOR 9
#define ADIOS_VERSION_PATCH 0

/* macros for comparing the version */
#define ADIOS_VERSION_GE(Maj,Min,Pat) \
       (((ADIOS_VERSION_MAJOR==Maj) && (ADIOS_VERSION_MINOR==Min) && (ADIOS_VERSION_PATCH>=Pat)) || \
        ((ADIOS_VERSION_MAJOR==Maj) && (ADIOS_VERSION_MINOR>Min)) || \
        (ADIOS_VERSION_MAJOR>Maj))

#define ADIOS_VERSION_LE(Maj,Min,Pat) \
       (((ADIOS_VERSION_MAJOR==Maj) && (ADIOS_VERSION_MINOR==Min) && (ADIOS_VERSION_PATCH<=Pat)) || \
        ((ADIOS_VERSION_MAJOR==Maj) && (ADIOS_VERSION_MINOR<Min)) || \
        (ADIOS_VERSION_MAJOR<Maj))

/* ADIOS Software release version as strings*/
#define ADIOS_VERSION_MAJOR_STRING "1"
#define ADIOS_VERSION_MINOR_STRING "9"
#define ADIOS_VERSION_PATCH_STRING "0"

/* BP File format version
 * BP versions 
 *   1 up until 1.6 release
 *   2 from 1.7 release: 32 bit variable IDs
 *   3 from 1.9 release: array attributes
 */
#define ADIOS_VERSION_BP_FORMAT 3

#endif
