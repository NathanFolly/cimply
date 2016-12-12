PETSC_ARCH = arch-linux-gnu-intel
CFLAGS = 
FFLAGS =
CPPFLAGS =
FPPFLAGS =
MANSEC =
LOCDIR =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

fimply: cimply.o fimply.o sim05reader.o 
	-${FLINKER} -o fimply cimply.o fimply.o sim05reader.o ${PETSC_LIB}
	${RM} cimply.o fimply.o sim05reader.o 

cimply: cimply.o chkopts
	-${CLINKER} -o cimply cimply.o ${PETSC_LIB}
	${RM} cimply.o

phcelltest: phcelltest.o phantomcell.o cimplyobjects.o
	-${CLINKER} -o phcelltest phcelltest.o phantomcell.o cimplyobjects.o ${PETSC_LIB}
	${RM} phcelltest.o phantomcell.o cimplyobjects.o

ooctest: ooctest.o phantomcell.o cimplyobjects.o testobject.o
	-${CLINKER} -o ooctest ooctest.o phantomcell.o cimplyobjects.o testobject.o ${PETSC_LIB}
	${RM} ooctest.o phantomcell.o cimplyobjects.o testobject.o
