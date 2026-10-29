# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Makefile for Rooster serial examples
#
# This file can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

# ATTENTION: path to SUNDIALS installation directory
prefix     = /usr/local
includedir = ${prefix}/fortran
libdir     = ${prefix}/lib

# ATTENTION: path to FORTRAN compiler
F90      = /usr/bin/gfortran
F90FLAGS = -O3 -DNDEBUG -O3 -DSUNDIALS_INT64_T -cpp
F90LIBS  =  -lm

# ATTENTION: build target
TARGET   = ZZBL
# ------------------------------------------------------------------------------

INCLUDES  = -I${includedir}
LIBRARIES = -lsundials_fida_mod -lsundials_ida -lsundials_fcore_mod -lsundials_core ${F90LIBS}
LINKFLAGS = -Wl,-rpath,$(libdir)

# ------------------------------------------------------------------------------

FILES    = ZZBL_fmod
FILES_DEPENDENCIES = C2_glodata_mod \
           C1_interpo_mod C1_readdat_mod \
           C0_corrsub_mod \
           A1_rescomp_mod \

OBJECTS =  ${FILES:=.o}
OBJECTS_DEPENDENCIES = ${FILES_DEPENDENCIES:=.o}

# ------------------------------------------------------------------------------

.SUFFIXES : .o .f90
.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $<

# ------------------------------------------------------------------------------

all: ${OBJECTS_DEPENDENCIES} ${OBJECTS}
	@for i in ${FILES} ; do \
	  echo "${F90} -o ${TARGET} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F90} -o ${TARGET} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

clean:
	rm -f *.o *.mod
	rm -f ${OBJECTS}
	rm -f ${OBJECTS_DEPENDENCIES}
	rm -f ${FILES}
	rm -f ${FILES_DEPENDENCIES}
# ------------------------------------------------------------------------------
