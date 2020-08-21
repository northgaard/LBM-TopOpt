#PETSC_DIR=/home/sebnorg/petsc-3.5.2
PETSC_DIR=/home/sebnorg/petsc-3.7.3
#PETSC_DIR=/home/sebnorg/petsc-3.5.4
#PETSC_ARCH=gcc-debug
PETSC_ARCH=gcc
LIBPATH=/home/sebnorg/Dropbox/PhD_MEK/Code/TOpTEnLB/src
#LIBPATH=/media/sebnorg/DataP1/Dropbox/PhD_MEK/Code/TOpTEnLB/src
CODIPATH=/home/sebnorg/CoDiPack-master/include
PETSCVER=-D PETSC35
TOPLIBARG=-lTOpTEnLB
#EXTCXX=-std=c++14
EXTCXX=-std=c++11
CXXFLAGS=-I. -I${LIBPATH} -I${CODIPATH} ${PETSCVER} ${EXTCXX}
FFLAGS=
FPPFLAGS=
LOCDIR=
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=
OBJFILES=diffuser.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

Diffuser2d: ${OBJFILES} chkopts
		rm -rf ./Diffuser2d
		-${CLINKER} -o ./Diffuser2d ${OBJFILES} ${PETSC_SYS_LIB} ${TOPLIBARG}
		${RM} ${OBJFILES}
		rm -rf *.o
