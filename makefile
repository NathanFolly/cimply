
CFLAGS = 
FFLAGS =
CPPFLAGS =
FPPFLAGS =
MANSEC =
LOCDIR =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

cimply: cimply.o chkopts
	-${CLINKER} -o cimply cimply.o ${PETSC_LIB}
	${RM} cimply.o
