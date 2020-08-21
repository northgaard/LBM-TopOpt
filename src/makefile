#PETSC_DIR=/home/sebnorg/petsc-3.5.2
#PETSC_DIR=/home/sebnorg/petsc-3.5.4
PETSC_DIR=/home/sebnorg/petsc-3.7.3
#PETSC_ARCH=gcc-debug
PETSC_ARCH=gcc
FULLPATH=${PETSC_DIR}/${PETSC_ARCH}/lib
LIBPATH=/home/sebnorg/Dropbox/PhD_MEK/Code/TOpTEnLB/src
CODIPATH=/home/sebnorg/CoDiPack-master/include
BOOSTPATH=/home/sebnorg/boost_1_62_0
PETSCVER=-D PETSC35
EXTCXX=-std=c++11
EXTC=-D GCCCOMP
CFLAGS=-I. -I${LIBPATH} ${PETSCVER} ${EXTC}
CPPFLAGS=-I. -I${LIBPATH} -I${CODIPATH} -I${BOOSTPATH} ${PETSCVER} ${EXTCXX}
FFLAGS=
FPPFLAGS=
LOCDIR=
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=
TOPLIBNAME=libTOpTEnLB
EXTERNFILES=external/Tapenade/adStack.o \
external/Tapenade/adBuffer.o \
external/SOBOL/sobol.o
OBJFILES=core/geometry2d.o \
core/memoryHandling.o \
core/MPISplit.o \
core/binaryIO.o \
latticeBoltzmann/lattices2d.o \
latticeBoltzmann/LBSolver.o \
latticeBoltzmann/NewLBSolver2d.o \
boundaryConditions/LBBoundaryLoop.o \
boundaryConditions/boundaryDescriptor2d.o \
boundaryConditions/bounceback2d.o \
adjointLatticeBoltzmann/adjointLBSolver.o \
adjointBoundaryConditions/adjointBounceback2d.o \
adjointBoundaryConditions/adjointSymmetry2d.o \
adjointBoundaryConditions/adjointZouHe2d.o \
topOpt/naiveCheckpointing.o \
topOpt/WangCheckpointing.o \
topOpt/staticCheckpointing.o \
topOpt/dynamicCheckpointing.o \
topOpt/LBCheckpointing.o \
topOpt/densityFilter2d.o \
topOpt/projectionFilter2d.o \
topOpt/variableProjectionFilter2d.o \
topOpt/doubleProjectionFilter2d.o \
topOpt/variableDoubleProjectionFilter2d.o \
topOpt/finiteDifferenceCheck2d.o \
topOpt/fullToDesignMapping2d.o \
topOpt/MMA.o \
functionals/volumeConstraint.o \
external/SobolGenerator.o \
core/geometry3d.o \
boundaryConditions/boundaryDescriptor3d.o \
latticeBoltzmann/LBSolver3d.o \
latticeBoltzmann/lattices3d.o
#latticeBoltzmann/multiGridSolver2d.o
#adjointBoundaryConditions/adjointBounceback3d.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

TOpTEnLB: ${EXTERNFILES} ${OBJFILES} chkopts
	-${CLINKER} -shared -Wl,-soname,${TOPLIBNAME}.so.1 -o \
	${FULLPATH}/${TOPLIBNAME}.so.1.0 ${OBJFILES} ${EXTERNFILES} ${PETSC_SYS_LIB}
	ln -sf ${FULLPATH}/${TOPLIBNAME}.so.1.0 ${FULLPATH}/${TOPLIBNAME}.so
	ln -sf ${FULLPATH}/${TOPLIBNAME}.so.1.0 ${FULLPATH}/${TOPLIBNAME}.so.1

myclean:
	rm ${OBJFILES}

allclean:
	rm ${OBJFILES} ${EXTERNFILES}
