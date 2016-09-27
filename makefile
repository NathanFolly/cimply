
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


