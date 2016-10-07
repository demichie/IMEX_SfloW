!*********************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!*********************************************************************
MODULE parameters

  IMPLICIT NONE

  REAL*8 :: eps_newton        !< threshold for the convergence of the
                              !< Newton's method 
  REAL*8 :: max_dt            !< Largest time step allowed
  REAL*8 :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL*8 :: eps_sing               !< parameter for desingularization

  REAL*8 :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

  !> Flag to add the relaxation terms after the linear reconstruction:\n
  !> - T      => evaluate the relaxation terms
  !> - F      => reconstruction without the relaxation 
  !> .
  LOGICAL :: interfaces_relaxation

  INTEGER, PARAMETER :: n_vars = 2        !< Number of conservative variables
  INTEGER, PARAMETER :: n_eqns = n_vars   !< Number of equations

  INTEGER :: n_nh     !< Number of non-hyperbolic terms

  INTEGER :: n_RK     !< Runge-Kutta order
  
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL*8, PARAMETER :: tol_abs = 1.D-5
  REAL*8, PARAMETER :: tol_rel = 1.D-5

  !> Limiter for the slope in the linear reconstruction:\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod sloe;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  INTEGER :: limiter(n_vars)

  !> Finite volume method:\n
  !> - 'LxF'       => lax-friedrichs scheme;
  !> - 'GFORCE '   => gforce scheme;
  !> - 'KT'        => Kurganov and Tadmor semidiscrete scheme;
  !> .
  CHARACTER(LEN=20) :: solver_scheme     

  REAL*8 :: theta             !< Van Leer limiter parameter
  REAL*8 :: t_start           !< initial time for the run
  REAL*8 :: t_end             !< end time for the run
  REAL*8 :: t_output          !< time of the next output
  REAL*8 :: dt_output         !< time interval for the output of the solution

  INTEGER :: verbose_level

  TYPE bc
     INTEGER :: flag
     REAL*8 :: value
  END TYPE bc

  ! -------boundary conditions variables

  !> blL&flag defines the left boundary condition:\n
  !> - bcL%flag = 0     => Dirichlet boundary condition;
  !> - bcL%flag = 1     => Neumann boundary condition.
  !> .
  !> bcL%value is the value of the left boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcL%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcL%flag=1).
  !> .
  TYPE(bc) :: bcL(n_vars)

  !> blR&flag defines the right boundary condition:\n
  !> - bcR%flag = 0     => Dirichlet boundary condition;
  !> - bcR%flag = 1     => Neumann boundary condition.
  !> .
  !> bcR%value is the value of the right boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcR%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcR%flag=1).
  !> .
  TYPE(bc) :: bcR(n_vars)

END MODULE parameters
