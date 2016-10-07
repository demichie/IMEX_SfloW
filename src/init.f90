!*********************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!*********************************************************************

MODULE init

  IMPLICIT none

  REAL*8, ALLOCATABLE :: q_init(:,:)
  
  !> Riemann problem interface relative position. It is a value
  !> between 0 and 1
  REAL*8 :: riemann_interface  


  REAL*8 :: hB_L          !< Left height
  REAL*8 :: u_L          !< Left velocity  
  
  REAL*8 :: hB_R          !< Right height
  REAL*8 :: u_R          !< Right velocity


CONTAINS


  !*********************************************************************
  !> \brief Riemann problem initialization
  !
  !> This subroutine initialize the solution for a Riemann problem. The 
  !> values for the left and right states and the interface location 
  !> are read from the input file.\
  !> \date 26/08/2011
  !*********************************************************************

  SUBROUTINE riemann_problem

    USE constitutive, ONLY : qp_to_qc

    USE geometry, ONLY : x0 , xN , x_comp , comp_cells , B_cent

    USE parameters, ONLY : n_vars , verbose_level

    USE solver, ONLY : q

    IMPLICIT none

    REAL*8 :: hB            !< height + batimetry
    REAL*8 :: u             !< velocity

    REAL*8 :: qp(n_vars) , qj(n_vars)

    INTEGER :: i            !< loop counter
    INTEGER :: i1           !< last index with left state

    REAL*8 :: eps

    riemann_int_search:DO i = 1,comp_cells

       IF ( x_comp(i) .LT. riemann_interface ) THEN

          i1 = i

       ELSE

          EXIT riemann_int_search

       END IF

    END DO riemann_int_search

    eps = 1.D-10

    ! Left initial state

    qp(1) = hB_L
    qp(2) = u_L


    DO i = 1,i1

       ! evaluate the vector of conservative variables
       CALL qp_to_qc( qp , B_cent(i) , qj )

       q(1:n_vars,i) = qj

       IF ( verbose_level .GE. 1 ) WRITE(*,*) i,B_cent(i),qp

    END DO

    ! Right initial state

    qp(1) = hB_R
    qp(2) = u_R


    DO i = i1+1,comp_cells

       ! evaluate the vector of conservative variables
       CALL qp_to_qc( qp , B_cent(i) , qj )

       q(1:n_vars,i) = qj

       IF ( verbose_level .GE. 1 ) WRITE(*,*) i,B_cent(i),qp
    
    END DO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE riemann_problem


END MODULE init
