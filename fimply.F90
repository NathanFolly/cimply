program fimply
  implicit none
#include "petsc/finclude/petscsys.h"
!  use iso_c_binding
  interface
     subroutine callthis() bind(c)
       use iso_c_binding
!       integer (c_int) :: argc
!       character (c_char) :: argv
     end subroutine callthis
  end interface
       
!integer (c_int) :: argc
!character (c_char) :: argv*1000
PetscErrorCode :: ierr

!argc = 1
!argv = "neumann"

!Since we cannot pass command line arguments in fortran we have to call the fortran version of the petscinitialize command that can cope without them. Options still can be set using an input file.
call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
call callthis()

call PetscFinalize(ierr);CHKERRQ(ierr)

end program fimply
