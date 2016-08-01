
CFLAGS = 
FFLAGS =
CPPFLAGS =
FPPFLAGS =
MANSEC =
LOCDIR =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

cimpleFEM: cimpleFEM.o chkopts
	-${CLINKER} -o cimpleFEM cimpleFEM.o ${PETSC_LIB}
	${RM} cimpleFEM.o
