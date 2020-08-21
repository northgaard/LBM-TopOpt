#PETSC_DIR=/home/sebnorg/petsc-3.5.2
PETSC_DIR=/home/sebnorg/petsc-3.7.3
#PETSC_DIR=/home/sebnorg/petsc-3.5.4
#PETSC_ARCH=gcc-debug
PETSC_ARCH=gcc
LIBPATH=/home/sebnorg/Dropbox/PhD_MEK/Code/TOpTEnLB/src
CODIPATH=/home/sebnorg/CoDiPack-master/include
PETSCVER=-D PETSC35
TOPLIBARG=-lTOpTEnLB
EXTCXX=-D GCCCOMP -std=c++11
CXXFLAGS=-I. -I${LIBPATH} -I${CODIPATH} ${PETSCVER} ${EXTCXX}
FFLAGS=
FPPFLAGS=
LOCDIR=
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=
OBJFILES=singleBlowHeatedObstacle.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

SingleBlowHeatedObstacle2d: ${OBJFILES} chkopts
		rm -rf ./SingleBlowHeatedObstacle2d
		-${CLINKER} -o ./SingleBlowHeatedObstacle2d ${OBJFILES} ${PETSC_SYS_LIB} ${TOPLIBARG}
		${RM} ${OBJFILES}
		rm -rf *.o
