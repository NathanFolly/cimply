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

interfacetest: interfacetest.o cimplyobjects.o interface.o simmsh.o femanalysis.o simmeranalysis.o
	-${CLINKER} -o interfacetest interfacetest.o cimplyobjects.o interface.o simmsh.o femanalysis.o simmeranalysis.o ${PETSC_LIB}
	${RM} cimplyobjects.o interface.o simmsh.o femanalysis.o simmeranalyis.o

boundaryvertextest: boundaryvertextest.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o
	-${CLINKER} -o boundaryvertextest boundaryvertextest.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o ${PETSC_LIB}
	${RM} boundaryvertextest.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o

phantommeshtest: phantommeshtest.o phantommesh.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o
	-${CLINKER} -o phantommeshtest phantommeshtest.o phantommesh.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o ${PETSC_LIB}
	${RM} phantommeshtest.o phantommesh.o cimplyobjects.o boundaryvertex.o geometry.o phantomcell.o simmsh.o interface.o simmeranalysis.o femanalysis.o
