!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!*********************************************************************
MODULE geometry

  USE parameters, ONLY : verbose_level

  IMPLICIT NONE

  !> Location of the centers of the control volume of the domain
  REAL*8, ALLOCATABLE :: x_comp(:)

  !> Location of the boundaries of the control volumes of the domain
  REAL*8, ALLOCATABLE :: x_stag(:)

  !> Batimetry at the boundaries of the control volumes
  REAL*8, ALLOCATABLE :: B_stag(:)

  !> Batimetry at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_cent(:)

  !> Batimetry slope at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_prime(:)

  REAL*8, ALLOCATABLE :: batimetry_profile(:,:)

  INTEGER :: n_batimetry_profile

  REAL*8 :: dx           !< Control volumes size
  REAL*8 :: x0           !< Left (bottom) of the physical domain
  REAL*8 :: xN           !< Right (top) of the physical domain
  INTEGER :: comp_cells  !< Number of control volumes in the comp. domain
  INTEGER :: comp_interfaces !< Number of interfaces (comp_cells+1)

CONTAINS

  !*********************************************************************
  !> \brief Finite volume grid initialization
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !*********************************************************************

  SUBROUTINE init_grid

    USE parameters, ONLY: eps_sing, batimetry_function_flag

    IMPLICIT none

    INTEGER j      !> loop counter

    comp_interfaces = comp_cells+1

    ALLOCATE( x_comp(comp_cells) )
    ALLOCATE( x_stag(comp_interfaces) )

    ALLOCATE( B_cent(comp_cells) )
    ALLOCATE( B_prime(comp_cells) )
    ALLOCATE( B_stag(comp_interfaces) )
    
    dx = ( xN - x0 ) / comp_cells
    
    eps_sing = MIN( dx , 1.D-5 ) ** 4

    WRITE(*,*) 'eps_sing = ',eps_sing

    x_comp(1) = x0 + 0.5D0 * dx
    x_stag(1) = x0

    ! rescaling in the case of batimetry defined in input file
    IF(batimetry_function_flag.EQV..FALSE.)THEN

       batimetry_profile(1,:) = x0 + ( xN - x0 ) * batimetry_profile(1,:) 

    ENDIF

    ! definition of batimetry by input file
    IF(batimetry_function_flag.EQV..FALSE.)THEN

      B_stag(1) = batimetry_profile(2,1)

      DO j=1,comp_cells

         x_stag(j+1) = x_stag(j) + dx
       
         x_comp(j) = 0.5 * ( x_stag(j) + x_stag(j+1) )

         CALL interp_1d_scalar( batimetry_profile(1,:) , batimetry_profile(2,:) , &
              x_stag(j+1) , B_stag(j+1) )


         B_cent(j) = 0.5 * ( B_stag(j) + B_stag(j+1) )
       
         B_prime(j) = ( B_stag(j+1) - B_stag(j) ) / (  x_stag(j+1) - x_stag(j) )

         IF ( verbose_level .GE. 2 ) THEN

            WRITE(*,*) batimetry_profile(1,:) 
            WRITE(*,*) batimetry_profile(2,:)
            WRITE(*,*) x_stag(j+1) , B_stag(j+1) ,x_comp(j) , B_cent(j) , B_prime(j) 
            READ(*,*)

         END IF

      END DO

    ! definition of batimetry by a function
    ELSE

      B_stag(1) = batimetry_function(x_stag(1))

      DO j=1,comp_cells

         x_stag(j+1) = x_stag(j) + dx
       
         x_comp(j) = 0.5 * ( x_stag(j) + x_stag(j+1) )

         B_stag(j+1) = batimetry_function(x_stag(j+1))

         B_cent(j) = 0.5 * ( B_stag(j) + B_stag(j+1) )
       
         B_prime(j) = ( B_stag(j+1) - B_stag(j) ) / (  x_stag(j+1) - x_stag(j) )

         IF ( verbose_level .GE. 2 ) THEN

            WRITE(*,*) batimetry_profile(1,:) 
            WRITE(*,*) batimetry_profile(2,:)
            WRITE(*,*) x_stag(j+1) , B_stag(j+1) ,x_comp(j) , B_cent(j) , B_prime(j) 
            READ(*,*)

         END IF

      END DO

    ENDIF

  END SUBROUTINE init_grid

!---------------------------------------------------------------------------
!> Scalar interpolation
!
!> This subroutine interpolate the values of the  array f1, defined on the 
!> grid points x1, at the point x2. The value are saved in f2
!> \date 13/02/2009
!> \param    x1           original grid                (\b input)
!> \param    f1           original values              (\b input)
!> \param    x2           new point                    (\b output)
!> \param    f2           interpolated value           (\b output)
!---------------------------------------------------------------------------

  SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
    IMPLICIT NONE
    
    REAL*8, INTENT(IN), DIMENSION(:) :: x1, f1
    REAL*8, INTENT(IN) :: x2
    REAL*8, INTENT(OUT) :: f2
    INTEGER :: n, n1x, t
    REAL*8 :: grad , rel_pos
    
    n1x = SIZE(x1)
  
    !
    ! ... locate the grid points near the topographic points
    ! ... and interpolate linearly the profile
    !
    t = 1

    search:DO n = 1, n1x-1

       rel_pos = ( x2 - x1(n) ) / ( x1(n+1) - x1(n) )

       IF ( ( rel_pos .GE. 0.D0 ) .AND. ( rel_pos .LE. 1.D0 ) ) THEN

          grad = ( f1(n+1)-f1(n) ) / ( x1(n+1)-x1(n) )
          f2 = f1(n) + ( x2-x1(n) ) * grad

          EXIT search
          
       ELSEIF  ( rel_pos .LT. 0.D0 ) THEN

          f2 = f1(n)
          
       ELSE

          f2 = f1(n+1)

       END IF

    END DO search

    RETURN
    
  END SUBROUTINE interp_1d_scalar


!---------------------------------------------------------------------------
!> Batimetry function
!
!> This subroutine generates a point of the batimetry from the input x grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!---------------------------------------------------------------------------
  REAL*8 FUNCTION batimetry_function(x)
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x

    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8, PARAMETER :: eps_dis = 10.0**(-8)

    ! example from Kurganov and Petrova 2007    
    IF(x.LT.0.0)THEN

      batimetry_function = 1.d0

    ELSEIF(x.GE.0.0.AND.x.LE.0.4)THEN

      batimetry_function = COS(pig*x)**2

    ELSEIF(x.GT.0.4.AND.x.LE.0.5)THEN

      batimetry_function = COS(pig*x)**2+0.25*(COS(10.0*pig*(x-0.5))+1)

    ELSEIF(x.GT.0.5.AND.x.LE.0.6)THEN

      batimetry_function = 0.5*COS(pig*x)**4+0.25*(COS(10.0*pig*(x-0.5))+1)

    ELSEIF(x.GT.0.6.AND.x.LT.1.0-eps_dis)THEN

      batimetry_function = 0.5*COS(pig*x)**4

    ELSEIF(x.GE.1.0-eps_dis.AND.x.LE.1.0+eps_dis)THEN

      batimetry_function = 0.25

    ELSEIF(x.GT.1.0+eps_dis.AND.x.LE.1.5)THEN

      batimetry_function = 0.25*SIN(2*pig*(x-1))

    ELSE

      batimetry_function = 0.d0

    ENDIF

  END FUNCTION batimetry_function
  
  
END MODULE geometry
