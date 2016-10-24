program fimply
  use iso_c_binding
  implicit none

#include "petsc/finclude/petscsys.h"

! Parameters
  integer, parameter :: NCALLS = 3  ! Number of times FEM is called
  integer, parameter :: IB=1         ! [SIMMER] N_cells in radial direction
  integer, parameter :: JB=11        ! [SIMMER] N_cells in vertical direction
  integer, parameter :: IBP2=IB+2    ! [SIMMER] IB plus shadow cells
  integer, parameter :: JBP2=JB+2    ! [SIMMER] JB plus shadow cells
  integer, parameter :: MMS = JBP2*IBP2! [SIMMER] total size of a scalar k-cell
                                       ! data array

  interface
     subroutine sim05tocimply() bind(c)
       use, intrinsic :: iso_c_binding
     end subroutine sim05tocimply
  end interface

  interface
     subroutine CimplySolve(MMS,PK,iter,tstep, ierr) bind(c)
       use,intrinsic :: iso_c_binding
       integer(c_int), value :: MMS
       PetscReal, Dimension(MMS) :: PK
       integer(c_int), value :: iter
       PetscReal, value :: tstep
      PetscErrorCode :: ierr
     end subroutine CimplySolve
  end interface

 interface
    subroutine CimplySetup(ierr) bind(c)
      use,intrinsic :: iso_c_binding
      PetscErrorCode :: ierr
    end subroutine CimplySetup
 end interface

 interface
    subroutine CimplyFinalize(ierr) bind(c)
      use, intrinsic :: iso_c_binding
      PetscErrorCode :: ierr
    end subroutine CimplyFinalize
 end interface

!Variables
integer :: m=0,n=0                       !counter
integer :: I=0                         !Radial position
integer :: J=0                         !vertical position
PetscReal :: tstep=5E-6
PetscReal, target :: PK(MMS)  ! The k-cell pressure

PetscErrorCode :: ierr

!Since we cannot pass command line arguments in fortran we have to call the fortran version of the petscinitialize command that can cope without them. Options still can be set using an input file.
call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
!call sim05tocimply()




!prepare the analysis

call CimplySetup(ierr);CHKERRQ(ierr);

! In SIMMER we call the FEM routine multiple times. 
! To try out how this works we use a do loop in this dummy code.

do m=1,NCALLS,1
   ! Assemble some kind of dummy k-cell pressure array
   do n=1,MMS ,1
      PK(n) = 0.0               !zero the pk array
   end do
   I=1
   do J=1,JB,1
      PK(IBP2*J+I+1)= 1E6*m+5E5*J
   end do
   call CimplySolve(MMS,PK,m,tstep,ierr);CHKERRQ(ierr);
end do


!clean up memory after the analysis

call CimplyFinalize(ierr);CHKERRQ(ierr);

call PetscFinalize(ierr);CHKERRQ(ierr)

end program fimply
