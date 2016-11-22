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
  REAL*8 :: u_L           !< Left velocity  
  
  REAL*8 :: hB_R          !< Right height
  REAL*8 :: u_R           !< Right velocity


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

    USE parameters, ONLY : n_vars , verbose_level, batimetry_function_flag

    USE solver, ONLY : q

    IMPLICIT none

    REAL*8 :: hB            !< height + batimetry
    REAL*8 :: u             !< velocity

    !REAL*8 :: qp(n_vars) , qj(n_vars)
    REAL*8 :: qp(n_vars,comp_cells) , qj(n_vars)

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

    DO i=1,i1

       qp(1,i) = hB_L

       qp(2,i) = u_L

    ENDDO

    DO i = 1,i1

       ! evaluate the vector of conservative variables
       !CALL qp_to_qc( qp , B_cent(i) , qj )
       CALL qp_to_qc( qp(:,i) , B_cent(i) , qj )

       q(1:n_vars,i) = qj

       IF ( verbose_level .GE. 1 ) WRITE(*,*) i,B_cent(i),qp

    END DO

    ! Right initial state

    DO i=i1+1,comp_cells

      qp(1,i) = hB_R

      qp(2,i) = u_R

    ENDDO

    DO i = i1+1,comp_cells

       ! evaluate the vector of conservative variables
       !CALL qp_to_qc( qp , B_cent(i) , qj )
       CALL qp_to_qc( qp(:,i) , B_cent(i) , qj )

       q(1:n_vars,i) = qj

       IF ( verbose_level .GE. 1 ) WRITE(*,*) i,B_cent(i),qp
    
    END DO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE riemann_problem

  !*********************************************************************
  !> \brief Problem initialization
  !
  !> This subroutine initialize the solution for a a generic problem. The 
  !> values are read from water and velocity functions
  !> \date OCTOBER 2016
  !*********************************************************************

  SUBROUTINE initial_conditions

    USE constitutive, ONLY : qp_to_qc

    USE geometry, ONLY : x0 , xN , x_comp , comp_cells , B_cent

    USE parameters, ONLY : n_vars , verbose_level, batimetry_function_flag

    USE solver, ONLY : q

    IMPLICIT none

    REAL*8 :: qp(n_vars,comp_cells) , qj(n_vars)

    INTEGER :: j          !< loop counter

    REAL*8 :: eps

    DO j=1,comp_cells

         qp(1,j) = water_function(x_comp(j),B_cent(j))

         qp(2,j) = velocity_function(x_comp(j),B_cent(j))

    ENDDO

    DO j = 1,comp_cells

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j) , B_cent(j) , qj )

         q(1:n_vars,j) = qj

         IF ( verbose_level .GE. 1 ) WRITE(*,*) j,B_cent(j),qp(:,j)

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE initial_conditions


!---------------------------------------------------------------------------
!> Water function
!
!> This subroutine defines water height hB=h+B
!> in the input x grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!---------------------------------------------------------------------------
  REAL*8 FUNCTION water_function(x,Bj)
    USE geometry, ONLY : batimetry_function

    IMPLICIT NONE
   
    REAL*8, INTENT(IN) :: x
    REAL*8, INTENT(IN) :: Bj

    ! example from Kurganov and Petrova 2007   
    !water_function=batimetry_function(x)
    water_function=Bj

    IF(x.LE.0.0)THEN
       water_function = water_function + 0.4d0
    ENDIF

  END FUNCTION water_function

!---------------------------------------------------------------------------
!> Velocity function
!
!> This subroutine defines the velocity
!> in the input x grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!---------------------------------------------------------------------------
  REAL*8 FUNCTION velocity_function(x,Bj)

    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x
    REAL*8, INTENT(IN) :: Bj
    
    ! example from Kurganov and Petrova 2007    
    velocity_function=0.d0

  END FUNCTION velocity_function


END MODULE init
