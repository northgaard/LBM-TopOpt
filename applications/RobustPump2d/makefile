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
#EXTCXX=-D GCCCOMP -std=c++14
EXTCXX=-std=c++1y
CXXFLAGS=-I. -I${LIBPATH} -I${CODIPATH} ${PETSCVER} ${EXTCXX}
FFLAGS=
FPPFLAGS=
LOCDIR=
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=
OBJFILES=robustPump2d.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

RobustPump2d: ${OBJFILES} chkopts
		rm -rf ./RobustPump2d
		-${CLINKER} -o ./RobustPump2d ${OBJFILES} ${PETSC_SYS_LIB} ${TOPLIBARG}
		${RM} ${OBJFILES}
		rm -rf *.o
