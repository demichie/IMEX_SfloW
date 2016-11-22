!*******************************************************************************
!> \brief Numerical solver
!
!> This module contains the variables and the subroutines for the 
!> numerical solution of the equations.  
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************
MODULE solver

  ! external variables

  USE constitutive, ONLY : implicit_flag

  USE geometry, ONLY : comp_cells
  USE geometry, ONLY : comp_interfaces

  USE geometry, ONLY : B_cent , B_prime , B_stag

  USE parameters, ONLY : n_eqns , n_vars , n_nh
  USE parameters, ONLY : n_RK
  USE parameters, ONLY : reconstr_variables
  USE parameters, ONLY : verbose_level

  USE parameters, ONLY : bcL , bcR

  IMPLICIT none
  
  !> Conservative variables
  REAL*8, ALLOCATABLE :: q(:,:)        
  !> Conservative variables
  REAL*8, ALLOCATABLE :: q0(:,:)        
  !> Solution of the finite-volume semidiscrete cheme
  REAL*8, ALLOCATABLE :: q_fv(:,:)     
  !> Reconstructed value at the left of the interface
  REAL*8, ALLOCATABLE :: q_interfaceL(:,:)        
  !> Reconstructed value at the right of the interface
  REAL*8, ALLOCATABLE :: q_interfaceR(:,:)
  !> Local speeds at the left of the interface
  REAL*8, ALLOCATABLE :: a_interfaceL(:,:)
  !> Local speeds at the right of the interface
  REAL*8, ALLOCATABLE :: a_interfaceR(:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL*8, ALLOCATABLE :: H_interface(:,:)
  !> Physical variables (\f$\alpha_1, p_1, p_2, \rho u, w, T\f$)
  REAL*8, ALLOCATABLE :: qp(:,:)
  

  !> Time step
  REAL*8 :: dt
  
  LOGICAL, ALLOCATABLE :: mask22(:,:) , mask21(:,:) , mask11(:,:) , mask12(:,:)

  !> Butcher Tableau for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_tilde_ij(:,:)
  !> Butcher Tableau for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_dirk_ij(:,:)

  !> Coefficients for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega_tilde(:)
  !> Coefficients for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega(:)

  !> Explicit coefficients for the hyperbolic part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_tilde(:)
  !> Explicit coefficients for the non-hyperbolic part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_dirk(:)
  !> Implicit coefficient for the non-hyperbolic part for a single step of the R-K scheme
  REAL*8 :: a_diag

  !> Intermediate solutions of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: q_rk(:,:,:)
  !> Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: F_x(:,:,:)
  !> Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NH(:,:,:)
 
  !> Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms(:,:,:)

  !> Local Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: Fxj(:,:)
  !> Local Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NHj(:,:)
  !> Local Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms_j(:,:)
  
  LOGICAL :: normalize_q
  LOGICAL :: normalize_f
  LOGICAL :: opt_search_NL

  REAL*8, ALLOCATABLE :: residual_term(:,:)


CONTAINS

  !*****************************************************************************
  !> \brief Memory allocation
  !
  !> This subroutine allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !*****************************************************************************

  SUBROUTINE allocate_solver_variables

    IMPLICIT NONE

    REAL*8 :: gamma

    INTEGER :: i,j

    ALLOCATE( q( n_vars , comp_cells ) , q0( n_vars , comp_cells ) )

    ALLOCATE( q_fv( n_vars , comp_cells ) )

    ALLOCATE( q_interfaceL( n_vars , comp_interfaces ) )
    ALLOCATE( q_interfaceR( n_vars , comp_interfaces ) )

    ALLOCATE( a_interfaceL( n_eqns , comp_interfaces ) )
    ALLOCATE( a_interfaceR( n_eqns , comp_interfaces ) )

    ALLOCATE( H_interface( n_eqns , comp_interfaces ) )

    ALLOCATE( qp( n_vars , comp_cells ) )

    ALLOCATE( a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( omega_tilde(n_RK) )
    ALLOCATE( omega(n_RK) )

    ALLOCATE( mask22(n_eqns,n_eqns) )
    ALLOCATE( mask21(n_eqns,n_eqns) )
    ALLOCATE( mask11(n_eqns,n_eqns) )
    ALLOCATE( mask12(n_eqns,n_eqns) )

    mask11(1:n_eqns,1:n_eqns) = .FALSE.
    mask12(1:n_eqns,1:n_eqns) = .FALSE.
    mask22(1:n_eqns,1:n_eqns) = .FALSE.
    mask21(1:n_eqns,1:n_eqns) = .FALSE.

    DO i = 1,n_eqns

       DO j = 1,n_eqns

          IF ( .NOT.implicit_flag(i) .AND. .NOT.implicit_flag(j) )              &
               mask11(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. .NOT.implicit_flag(j) )                   &
               mask12(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. implicit_flag(j) )                        &
               mask22(j,i) = .TRUE.
          IF ( .NOT.implicit_flag(i) .AND. implicit_flag(j) )                   &
               mask21(j,i) = .TRUE.

       END DO

    END DO


    a_tilde_ij = 0.D0
    a_dirk_ij = 0.D0
    omega_tilde = 0.D0
    omega = 0.D0

    gamma = 1.D0 - 1.D0 / SQRT(2.D0)

    IF ( n_RK .EQ. 1 ) THEN

       a_tilde_ij(1,1) = 1.D0
       
       omega_tilde(1) = 1.D0

       a_dirk_ij(1,1) = 0.D0

       omega(1) = 0.D0

    ELSEIF ( n_RK .EQ. 2 ) THEN

       a_tilde_ij(2,1) = 1.0D0

       omega_tilde(1) = 1.0D0
       omega_tilde(2) = 0.0D0

       a_dirk_ij(2,2) = 1.0D0
        
       omega(1) = 0.D0
       omega(2) = 1.D0

    ELSEIF ( n_RK .EQ. 3 ) THEN

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 0.5D0
       a_tilde_ij(3,2) = 0.5D0

       omega_tilde(1) =  1.0D0 / 3.0D0
       omega_tilde(2) =  1.0D0 / 3.0D0
       omega_tilde(3) =  1.0D0 / 3.0D0

       a_dirk_ij(1,1) = 0.25D0
       a_dirk_ij(2,2) = 0.25D0
       a_dirk_ij(3,1) = 1.0D0 / 3.0D0
       a_dirk_ij(3,2) = 1.0D0 / 3.0D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0
        
       omega(1) =  1.0D0 / 3.0D0
       omega(2) =  1.0D0 / 3.0D0
       omega(3) =  1.0D0 / 3.0D0

    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! LRR(3,2,2)

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 1.D0 / 3.D0
       a_tilde_ij(4,2) = 1.0D0

       omega_tilde(1) = 0.D0
       omega_tilde(2) = 1.0D0
       omega_tilde(3) = 0.0D0
       omega_tilde(4) = 0.D0

       a_dirk_ij(2,2) = 0.5D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0
       a_dirk_ij(4,3) = 0.75D0
       a_dirk_ij(4,4) = 0.25D0
        
       omega(1) = 0.D0
       omega(2) = 0.D0
       omega(3) = 0.75D0
       omega(4) = 0.25D0

    END IF

    ALLOCATE( a_tilde(n_RK) )
    ALLOCATE( a_dirk(n_RK) )
    
    ALLOCATE( q_rk( n_vars , comp_cells , n_RK ) )
    ALLOCATE( F_x( n_eqns , comp_cells , n_RK ) )
    ALLOCATE( NH( n_eqns , comp_cells , n_RK ) )

    ALLOCATE( expl_terms( n_eqns , comp_cells , n_RK ) )

    ALLOCATE( Fxj(n_eqns,n_RK) )
    ALLOCATE( NHj(n_eqns,n_RK) )
    ALLOCATE( expl_terms_j(n_eqns,n_RK) )

    ALLOCATE( residual_term( n_vars , comp_cells ) )

  END SUBROUTINE allocate_solver_variables

  !******************************************************************************
  !> \brief Memory deallocation
  !
  !> This subroutine de-allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE deallocate_solver_variables

    DEALLOCATE( q , q0 )

    DEALLOCATE( q_fv )

    DEALLOCATE( q_interfaceL )
    DEALLOCATE( q_interfaceR )

    DEALLOCATE( a_interfaceL )
    DEALLOCATE( a_interfaceR )

    DEALLOCATE( H_interface )

    Deallocate( qp )

    DEALLOCATE( a_tilde_ij )
    DEALLOCATE( a_dirk_ij )
    DEALLOCATE( omega_tilde )
    DEALLOCATE( omega )

    DEALLOCATE( implicit_flag )

    DEALLOCATE( a_tilde )
    DEALLOCATE( a_dirk )
   
    DEALLOCATE( q_rk )
    DEALLOCATE( F_x )
    DEALLOCATE( NH )

    DEALLOCATE( Fxj )
    DEALLOCATE( NHj )

    DEALLOCATE( mask22 , mask21 , mask11 , mask12 )

    DEALLOCATE( residual_term )

  END SUBROUTINE deallocate_solver_variables

  !*****************************************************************************
  !> \brief Time-step computation
  !
  !> This subroutine evaluate the maximum time step according to the CFL
  !> condition. The local speed are evaluated with the characteristic
  !> polynomial of the Jacobian of the fluxes.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !*****************************************************************************

  SUBROUTINE timestep

    ! External variables
    USE geometry, ONLY : dx
    USE parameters, ONLY : max_dt , cfl

    ! External procedures
    USE constitutive, ONLY : eval_local_speeds

    IMPLICIT none

    REAL*8 :: vel_max(n_vars)
    REAL*8 :: vel_min(n_vars)
    REAL*8 :: vel_j         !< maximum speed in the j-th cell
    REAL*8 :: dt_cfl        !< local time step
    REAL*8 :: qj(n_vars)    !< conservative variables

    REAL*8 :: a_interface_max(n_eqns,comp_interfaces)
    REAL*8 :: dt_interface

    INTEGER :: i, j            !< loop counters

    dt = max_dt

    IF ( cfl .NE. -1.d0 ) THEN

       IF ( reconstr_variables .EQ. 0 ) THEN

          ! Linear reconstruction of the conservative variables at the interfaces
          CALL reconstruction_cons

       ELSEIF ( reconstr_variables .EQ. 0 ) THEN
     
          ! Linear reconstruction of the physical var. (h+B,u) at the interfaces
          CALL reconstruction_phys1
        
       ELSE 

          ! Linear reconstruction of the physical var. (h,u) at the interfaces
          CALL reconstruction_phys2
      
       END IF

       ! Evaluation of the maximum local speeds at the interfaces
       CALL eval_speeds
  
       DO i=1,n_vars

          a_interface_max(i,:) = MAX( a_interfaceR(i,:) , -a_interfaceL(i,:) )

       END DO

       DO j = 1,comp_cells

          ! qj = q( 1:n_vars , j )
          
          ! CALL eval_local_speeds(qj,B_cent(j),vel_min,vel_max)

          ! vel_j = MAX( MAXVAL(ABS(vel_min)) , MAXVAL(ABS(vel_max)) )
          
          ! dt_cfl = cfl * dx / vel_j
          
          ! dt = MIN( dt , dt_cfl )
     
          dt_interface = cfl * dx / MAX( MAXVAL(a_interface_max(:,j)) ,         &
               MAXVAL(a_interface_max(:,j+1)) )

          dt = MIN( dt , dt_interface )
     
       END DO
     
    END IF




  END SUBROUTINE timestep

  !*****************************************************************************
  !> \brief Runge-Kutta integration
  !
  !> This subroutine integrate the hyperbolic conservation law with
  !> non-hyperbolic terms using an implicit-explicit runge-kutta scheme.
  !> The fluxes are integrated explicitely while the non-hyperbolic terms
  !> are integrated implicitely.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************
  
  SUBROUTINE imex_RK_solver
    
    USE constitutive, ONLY : eval_nonhyperbolic_terms
    
    USE constitutive, ONLY : qc_to_qp
    
    IMPLICIT NONE
    
    REAL*8 :: q_guess(n_vars) !< initial guess for the solution of the RK step
    INTEGER :: i_RK           !< loop counter for the RK iteration
    INTEGER :: j              !< loop counter over the grid volumes
    
    REAL*8 :: nh_term(n_eqns) !< non_hyperbolic terms

    INTEGER :: i

    q0( 1:n_vars , 1:comp_cells )  = q( 1:n_vars , 1:comp_cells )

    IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: beginning'
    
    ! Initialization of the variables for the Runge-Kutta scheme
    q_rk(1:n_vars,1:comp_cells,1:n_RK) = 0.d0
    
    F_x(1:n_eqns,1:comp_cells,1:n_RK) = 0.d0
    
    NH(1:n_eqns,1:comp_cells,1:n_RK) = 0.d0
    
    expl_terms(1:n_eqns,1:comp_cells,1:n_RK) = 0.d0
    
    DO i_RK = 1,n_RK
       
       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       a_tilde = 0.d0
       a_dirk = 0.d0
       
       a_tilde(1:i_RK-1) = a_tilde_ij(i_RK,1:i_RK-1)
       a_dirk(1:i_RK-1) = a_dirk_ij(i_RK,1:i_RK-1)
       
       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = a_dirk_ij(i_RK,i_RK)

       DO j = 1,comp_cells

          IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: j',j
          
          IF ( i_RK .EQ. 1 ) THEN
             
             q_guess(1:n_vars) = q0( 1:n_vars , j ) 

          ELSE
             
             q_guess(1:n_vars) = q_rk( 1:n_vars , j , MAX(1,i_RK-1) )
             
          END IF

          Fxj(1:n_eqns,1:n_RK) = F_x( 1:n_eqns , j , 1:n_RK )

          NHj(1:n_eqns,1:n_RK) = NH( 1:n_eqns , j , 1:n_RK )

          Expl_terms_j(1:n_eqns,1:n_RK) = expl_terms( 1:n_eqns , j , 1:n_RK )

          IF ( verbose_level .GE. 2 ) THEN
             
             WRITE(*,*) 'q_guess',q_guess
             CALL qc_to_qp( q_guess , B_cent(j) , qp(1:n_vars,j) )
             WRITE(*,*) 'q_guess: qp',qp(1:n_vars,j)

          END IF

          IF ( a_diag .NE. 0.D0 ) THEN

             CALL solve_rk_step( B_cent(j) , B_prime(j) , q_guess ,             &
                  q0(1:n_vars,j ) , a_tilde , a_dirk , a_diag )

          END IF

          q_rk( 1:n_vars , j , i_RK ) = q_guess


          ! store the non-hyperbolic term for the explicit computations
          IF ( a_diag .EQ. 0.D0 ) THEN

             CALL eval_nonhyperbolic_terms( B_cent(j) , B_prime(j) ,            &
                  r_qj = q_guess , r_nh_term_impl = NH(1:n_eqns,j,i_RK) ) 

          ELSE
         
             NH( 1:n_eqns , j , i_RK ) = 1.D0 / a_diag * ( ( q_guess -          &
                  q0( 1:n_vars , j ) ) / dt +                                   &
                  ( MATMUL(Fxj,a_tilde) - MATMUL(NHj,a_dirk) ) )

          END IF

          IF ( verbose_level .GE. 2 ) THEN
             
             WRITE(*,*) 'imex_RK_solver: qc',q_guess
             CALL qc_to_qp( q_guess, B_cent(j) , qp(1:n_vars,j) )
             WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars,j)
             READ(*,*)

          END IF

       END DO

       ! Eval and save the explicit hyperbolic (fluxes) terms
       CALL eval_hyperbolic_terms( q_rk(1:n_vars,1:comp_cells,i_RK) ,           &
            F_x(1:n_eqns,1:comp_cells,i_RK) )
       
       ! Eval and save the other explicit terms (e.g. gravity or viscous forces)
       CALL eval_explicit_terms( q_rk(1:n_vars,:,i_RK) ,                        &
            expl_terms(1:n_eqns,:,i_RK) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'Fx, expl_terms'

          DO j = 1,comp_cells
             
             WRITE(*,*) F_x(1:n_vars,j,i_RK),expl_terms(2,j,i_RK)
             
          END DO
          
          READ(*,*)

       END IF

    END DO

    DO j = 1,comp_cells
       
       residual_term(1:n_vars,j) = MATMUL( F_x(1:n_eqns,j,1:n_RK)               &
            + expl_terms(1:n_eqns,j,1:n_RK) , omega_tilde )                     &
            - MATMUL( NH(1:n_eqns,j,1:n_RK) , omega )

    END DO
    


    DO j = 1,comp_cells

       IF ( verbose_level .GE. 1 ) THEN
          
          WRITE(*,*) 'cell j =',j
          WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j)
          CALL qc_to_qp(q0(1:n_vars,j) , B_cent(j) , qp(1:n_vars,j))
          WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j)

       END IF

       q(1:n_vars,j) = q0(1:n_vars,j) - dt * residual_term(1:n_vars,j)

       IF ( verbose_level .GE. 1 ) THEN

          CALL qc_to_qp(q(1:n_vars,j) , B_cent(j) , qp(1:n_vars,j))

          WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j)
          WRITE(*,*) 'after imex_RK_solver: qp',qp(1:n_vars,j)
          READ(*,*)

       END IF

    END DO

  END SUBROUTINE imex_RK_solver
  
  !******************************************************************************
  !> \brief Runge-Kutta single step integration
  !
  !> This subroutine find the solution of the non-linear system 
  !> given the a step of the implicit-explicit Runge-Kutta scheme for a
  !> cell:\n
  !> \f$ Q^{(i)} = Q^n - dt \sum_{j=1}^{i-1}\tilde{a}_{j}\partial_x 
  !> F(Q^{(j)}) +  dt \sum_{j=1}^{i-1} a_j  NH(Q^{(j)}) 
  !> + dt a_{diag} NH(Q^{(i)}) \f$\n
  !
  !> \param[in]     Bj        batimetry at the cell center
  !> \param[in]     Bprimej   batimetry slope at the cell center
  !> \param[in,out] qj        conservative variables 
  !> \param[in]     qj_old    conservative variables at the old time step
  !> \param[in]     a_tilde   explicit coefficents for the fluxes
  !> \param[in]     a_dirk    explicit coefficient for the non-hyperbolic terms
  !> \param[in]     a_diag    implicit coefficient for the non-hyperbolic terms 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************
  
  SUBROUTINE solve_rk_step( Bj , Bprimej, qj, qj_old, a_tilde , a_dirk , a_diag )

    USE parameters, ONLY : max_nl_iter , tol_rel , tol_abs

    USE constitutive, ONLY : qc_to_qp

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej
    REAL*8, INTENT(INOUT) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag

    REAL*8 :: qj_org(n_vars) , qj_rel(n_vars)
    
    REAL*8 :: left_matrix(n_eqns,n_vars)
    REAL*8 :: right_term(n_eqns)

    REAL*8 :: scal_f

    REAL*8 :: coeff_f(n_eqns)

    REAL*8 :: qj_rel_NR_old(n_vars)
    REAL*8 :: scal_f_old
    REAL*8 :: desc_dir(n_vars)
    REAL*8 :: grad_f(n_vars)

    INTEGER :: pivot(n_vars)

    REAL*8 :: desc_dir_small(n_vars)

    REAL*8 :: left_matrix_small22(n_nh,n_nh)
    REAL*8 :: left_matrix_small21(n_eqns-n_nh,n_nh)
    REAL*8 :: left_matrix_small11(n_eqns-n_nh,n_vars-n_nh)
    REAL*8 :: left_matrix_small12(n_nh,n_vars-n_nh)

    REAL*8 :: right_term_small2(n_nh)
    REAL*8 :: desc_dir_small2(n_nh)
    INTEGER :: pivot_small2(n_nh)

    REAL*8 :: right_term_small1(n_eqns-n_nh)
    REAL*8 :: desc_dir_small1(n_vars-n_nh)
    INTEGER :: pivot_small1(n_vars-n_nh)

    INTEGER :: ok
    
    INTEGER :: i 
    INTEGER :: nl_iter

    REAL*8, PARAMETER :: STPMX=100.D0
    REAL*8 :: stpmax
    LOGICAL :: check

    REAL*8, PARAMETER :: TOLF=1.D-10 , TOLMIN=1.D-6
    REAL*8 :: TOLX

    REAL*8 :: scal_f_init

    REAL*8 :: qpj(n_vars)

    REAL*8 :: desc_dir2(n_vars)

    REAL*8 :: desc_dir_temp(n_vars)

    normalize_q = .TRUE.
    normalize_f = .FALSE.
    opt_search_NL = .TRUE.
    
    coeff_f(1:n_eqns) = 1.D0

    ! normalize the functions of the nonlinear system

    IF ( normalize_f ) THEN
       
       qj = qj_old - dt * ( MATMUL(Fxj,a_tilde) - MATMUL(NHj,a_dirk) )

       CALL eval_f( Bj , Bprimej , qj , qj_old , a_tilde , a_dirk , a_diag ,    &
            coeff_f , right_term , scal_f )

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'solve_rk_step: non-normalized right_term'
          WRITE(*,*) right_term
          WRITE(*,*) 'scal_f',scal_f

       END IF
        
       DO i=1,n_eqns
          
          IF ( ABS(right_term(i)) .GE. 1.D0 ) coeff_f(i) = 1.D0 / right_term(i)

       END DO

       right_term = coeff_f * right_term
       
       scal_f = 0.5D0 * DOT_PRODUCT( right_term , right_term )

       IF ( verbose_level .GE. 3 ) THEN                    
          WRITE(*,*) 'solve_rk_step: after normalization',scal_f
       END IF

    END IF

    ! set the initial guess for the NR iterative solver

!!$    qj = qj_old - dt * ( MATMUL(Fxj,a_tilde) - MATMUL(NHj,a_dirk) )
!!$
!!$    DO i=1,n_eqns
!!$       
!!$       IF ( implicit_flag(i) ) qj(i) = qj_old(i)
!!$       
!!$    END DO

!    qj = qj_old

    !---- normalize the conservative variables ------

    IF ( normalize_q ) THEN

       qj_org = qj
       
       qj_org = MAX( ABS(qj_org) , 1.D-3 )

    ELSE 

       qj_org(1:n_vars) = 1.D0

    END IF
    
    qj_rel = qj / qj_org

    ! -----------------------------------------------

    DO nl_iter=1,max_nl_iter

       TOLX = epsilon(qj_rel)
       
       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solve_rk_step: nl_iter',nl_iter
       
       CALL eval_f( Bj , Bprimej , qj , qj_old , a_tilde , a_dirk , a_diag ,    &
            coeff_f , right_term , scal_f )
       
       IF ( verbose_level .GE. 2 ) THEN
          
          WRITE(*,*) 'solve_rk_step: right_term',right_term
       
       END IF

       IF ( verbose_level .GE. 2 ) THEN
          
          WRITE(*,*) 'before_lnsrch: scal_f',scal_f
       
       END IF

       ! check the residual of the system
       
       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN
          
          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          RETURN
          
       END IF
       
       IF ( ( normalize_f ) .AND. ( scal_f < 1.D-6 ) ) THEN
          
          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check scal_f',check
          RETURN
          
       END IF

       ! ---- evaluate the descent direction ------------------------------------
    
       IF ( COUNT( implicit_flag ) .EQ. n_eqns ) THEN

          CALL eval_jacobian(Bj,Bprimej,qj_rel,qj_org,coeff_f,left_matrix)
          
          desc_dir_temp = - right_term
          
          CALL DGESV(n_eqns,1, left_matrix , n_eqns, pivot, desc_dir_temp ,     &
               n_eqns, ok)
          
          desc_dir = desc_dir_temp

       ELSE

          CALL eval_jacobian(Bj,Bprimej,qj_rel,qj_org,coeff_f,left_matrix)
          
          left_matrix_small11 = reshape(pack(left_matrix, mask11),              &
               [n_eqns-n_nh,n_eqns-n_nh]) 
          
          left_matrix_small12 = reshape(pack(left_matrix, mask12),              &
               [n_nh,n_eqns-n_nh]) 
          
          left_matrix_small22 = reshape(pack(left_matrix, mask22),              &
               [n_nh,n_nh]) 
          
          left_matrix_small21 = reshape(pack(left_matrix, mask21),              &
               [n_eqns-n_nh,n_nh]) 
          
          
          desc_dir_small1 = pack( right_term, .NOT.implicit_flag )
          desc_dir_small2 = pack( right_term , implicit_flag )
          
          DO i=1,n_vars-n_nh
             
             desc_dir_small1(i) = desc_dir_small1(i) / left_matrix_small11(i,i)
             
          END DO
          
          desc_dir_small2 = desc_dir_small2 -                                   &
               MATMUL( desc_dir_small1 , left_matrix_small21 )
          
          CALL DGESV(n_nh,1, left_matrix_small22 , n_nh , pivot_small2 ,        &
               desc_dir_small2 , n_nh, ok)
          
          desc_dir = unpack( - desc_dir_small2 , implicit_flag , 0.0D0 )        &
               + unpack( - desc_dir_small1 , .NOT.implicit_flag , 0.0D0 )
          
       END IF

       IF ( verbose_level .GE. 3 ) WRITE(*,*) 'desc_dir',desc_dir

       qj_rel_NR_old = qj_rel
       scal_f_old = scal_f

       IF ( ( opt_search_NL ) .AND. ( nl_iter .GT. 1 ) ) THEN
       ! Search for the step lambda giving a sufficient decrease in the solution 

          stpmax = STPMX * MAX( SQRT( DOT_PRODUCT(qj_rel,qj_rel) ) ,            &
               DBLE( SIZE(qj_rel) ) )

          grad_f = MATMUL( right_term , left_matrix )

          desc_dir2 = desc_dir

          CALL lnsrch( Bj , Bprimej , qj_rel_NR_old , qj_org , qj_old ,         &
               scal_f_old , grad_f , desc_dir , coeff_f , qj_rel , scal_f ,     &
               right_term , stpmax , check )
          
       ELSE
          
          qj_rel = qj_rel_NR_old + desc_dir
          
          qj = qj_rel * qj_org
          
          CALL eval_f( Bj , Bprimej , qj , qj_old , a_tilde , a_dirk , a_diag , &
               coeff_f , right_term , scal_f )

       END IF

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'after_lnsrch: scal_f',scal_f

       qj = qj_rel * qj_org

       IF ( verbose_level .GE. 3 ) THEN
       
          WRITE(*,*) 'qj',qj
          CALL qc_to_qp( qj , Bj , qpj)
          WRITE(*,*) 'qp',qpj

       END IF

       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          check= .FALSE.
          RETURN

       END IF
		
       IF (check) THEN

          check = ( MAXVAL( ABS(grad_f(:)) * MAX( ABS( qj_rel(:) ),1.D0 ) /     &
               MAX( scal_f , 0.5D0 * SIZE(qj_rel) ) )  < TOLMIN )

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '2: check',check
!          RETURN
	
       END IF

       IF ( MAXVAL( ABS( qj_rel(:) - qj_rel_NR_old(:) ) / MAX( ABS( qj_rel(:)) ,&
            1.D0 ) ) < TOLX ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check',check
          RETURN

       END IF

    END DO

  END SUBROUTINE solve_rk_step

  !******************************************************************************
  !> \brief Search the descent stepsize
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> \param[in]     Bj             batimetry at the cell center
  !> \param[in]     Bprimej        batimetry slope at the cell center
  !> \param[in]     qj_rel_NR_old  
  !> \param[in]     qj_org
  !> \param[in]     qj_old
  !> \param[in]     scal_f_old
  !> \param[in]     grad_f
  !> \param[in,out] desc_dir
  !> \param[in]     coeff_f
  !> \param[out]    qj_rel
  !> \param[out]    scal_f
  !> \param[out]    right_term
  !> \param[in]     stpmax
  !> \param[out]    check
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE lnsrch( Bj , Bprimej , qj_rel_NR_old , qj_org , qj_old ,           &
       scal_f_old , grad_f , desc_dir , coeff_f , qj_rel , scal_f , right_term ,&
       stpmax , check )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj

    REAL*8, INTENT(IN) :: Bprimej

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_org

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_old

    !> Gradient at xold
    REAL*8, DIMENSION(:), INTENT(IN) :: grad_f

    !> Value of the function at xold
    REAL*8, INTENT(IN) :: scal_f_old

    !> Descent direction (usually Newton direction)
    REAL*8, DIMENSION(:), INTENT(INOUT) :: desc_dir

    REAL*8, INTENT(IN) :: stpmax

    !> Coefficients to rescale the nonlinear function
    REAL*8, DIMENSION(:), INTENT(IN) :: coeff_f

    !> Updated solution
    REAL*8, DIMENSION(:), INTENT(OUT) :: qj_rel

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: scal_f

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: right_term(n_eqns)

    !> Output quantity check is false on a normal exit 
    LOGICAL, INTENT(OUT) :: check

    REAL*8, PARAMETER :: TOLX=epsilon(qj_rel)

    INTEGER, DIMENSION(1) :: ndum
    REAL*8 :: ALF , a,alam,alam2,alamin,b,disc
    REAL*8 :: scal_f2
    REAL*8 :: desc_dir_abs
    REAL*8 :: rhs1 , rhs2 , slope, tmplam

    REAL*8 :: scal_f_min , alam_min

    REAL*8 :: qj(n_vars)

    ALF = 1.0d-4

    IF ( size(grad_f) == size(desc_dir) .AND. size(grad_f) == size(qj_rel) .AND.    &
         size(qj_rel) == size(qj_rel_NR_old) ) THEN

       ndum = size(grad_f)

    ELSE

       WRITE(*,*) 'nrerror: an assert_eq failed with this tag:', 'lnsrch'
       STOP 'program terminated by assert_eq4'

    END IF

    check = .FALSE.

    desc_dir_abs = SQRT( DOT_PRODUCT(desc_dir,desc_dir) )

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax / desc_dir_abs  
    
    slope = DOT_PRODUCT(grad_f,desc_dir)

    alamin = TOLX / MAXVAL( ABS( desc_dir(:)) / MAX( ABS(qj_rel_NR_old(:)),1.D0 ) )

    IF ( alamin .EQ. 0.d0) THEN

       qj_rel(:) = qj_rel_NR_old(:)

       RETURN

    END IF

    alam = 1.0D0
    
    scal_f_min = scal_f_old

    optimal_step_search: DO

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'alam',alam

       END IF

       qj_rel = qj_rel_NR_old + alam * desc_dir
       
       qj = qj_rel * qj_org

       CALL eval_f( Bj , Bprimej , qj , qj_old , a_tilde , a_dirk , a_diag ,    &
            coeff_f , right_term , scal_f )

       IF ( verbose_level .GE. 4 ) THEN
          
          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          READ(*,*)
     
       END IF

       IF ( scal_f .LT. scal_f_min ) THEN

          scal_f_min = scal_f
          alam_min = alam
          
       END IF

       IF ( scal_f .LE. 0.9 * scal_f_old ) THEN   
          ! sufficient function decrease
          
          IF ( verbose_level .GE. 4 ) THEN
             
             WRITE(*,*) 'sufficient function decrease'
             
          END IF
          
          EXIT optimal_step_search   

       ELSE IF ( alam < alamin ) THEN   
          ! convergence on Delta_x

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) ' convergence on Delta_x',alam,alamin

          END IF
          
          qj_rel(:) = qj_rel_NR_old(:)
          scal_f = scal_f_old
          check = .TRUE.
          
          EXIT optimal_step_search

!       ELSE IF ( scal_f .LE. scal_f_old + ALF * alam * slope ) THEN   
       ELSE  

          IF ( alam .EQ. 1.D0 ) THEN

             tmplam = - slope / ( 2.0D0 * ( scal_f - scal_f_old - slope ) )

          ELSE

             rhs1 = scal_f - scal_f_old - alam*slope
             rhs2 = scal_f2 - scal_f_old - alam2*slope

             a = ( rhs1/alam**2.D0 - rhs2/alam2**2.D0 ) / ( alam - alam2 )
             b = ( -alam2*rhs1/alam**2 + alam*rhs2/alam2**2 ) / ( alam - alam2 )
             
             IF ( a .EQ. 0.D0 ) THEN

                tmplam = - slope / ( 2.0D0 * b )

             ELSE

                disc = b*b - 3.0D0*a*slope

                IF ( disc .LT. 0.D0 ) THEN

                   tmplam = 0.5D0 * alam

                ELSE IF ( b .LE. 0.D0 ) THEN

                   tmplam = ( - b + SQRT(disc) ) / ( 3.D0 * a )

                ELSE

                   tmplam = - slope / ( b + SQRT(disc) )

                ENDIF

             END IF

             IF ( tmplam .GT. 0.5D0*alam ) tmplam = 0.5D0 * alam

          END IF

       END IF

       alam2 = alam
       scal_f2 = scal_f
       alam = MAX( tmplam , 0.5D0*alam )

    END DO optimal_step_search
    
  END SUBROUTINE lnsrch
  
  !******************************************************************************
  !> \brief Evaluate the nonlinear system
  !
  !> This subroutine evaluate the value of the nonlinear system in the state 
  !> defined by the variables qj.
  !> \param[in]    Bj          batimetry at the cell center
  !> \param[in]    Bprimej     batimetry slope at the cell center
  !> \param[in]    qj          conservative variables 
  !> \param[in]    qj_old      conservative variables at the old time step
  !> \param[in]    a_tilde     explicit coefficients for the hyperbolic terms 
  !> \param[in]    a_dirk      explicit coefficients for the non-hyperbolic terms 
  !> \param[in]    a_diag      implicit coefficient for the non-hyperbolic term
  !> \param[in]    coeff_f     coefficient to rescale the nonlinear functions
  !> \param[out]   f_nl        values of the nonlinear functions
  !> \param[out]   scal_f      value of the scalar function f=0.5*<F,F>
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_f( Bj , Bprimej , qj , qj_old , a_tilde , a_dirk , a_diag ,   &
       coeff_f , f_nl , scal_f )
  
    USE constitutive, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej
    REAL*8, INTENT(IN) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: f_nl(n_eqns)
    REAL*8, INTENT(OUT) :: scal_f

    REAL*8 :: nh_term_impl(n_eqns)
    REAL*8 :: Rj(n_eqns)

    CALL eval_nonhyperbolic_terms( Bj , Bprimej , r_qj = qj ,                   &
         r_nh_term_impl = nh_term_impl ) 

    Rj = ( MATMUL(Fxj,a_tilde) - MATMUL(NHj,a_dirk) ) - a_diag * nh_term_impl
    
    f_nl = qj - qj_old + dt * Rj
        
    f_nl = coeff_f * f_nl

    scal_f = 0.5D0 * DOT_PRODUCT( f_nl , f_nl )

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Evaluate the jacobian 
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !> \param[in]    Bj          batimetry at the cell center
  !> \param[in]    Bprimej     batimetry slope at the cell center
  !> \param[in]    qj_rel      relative variation (qj=qj_rel*qj_org)
  !> \param[in]    qj_org      conservative variables at the old time step
  !> \param[in]    coeff_f     coefficient to rescale the nonlinear functions
  !> \param[out]   left_matrix matrix from the linearization of the system
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_jacobian(Bj , Bprimej , qj_rel , qj_org , coeff_f, left_matrix)

    USE constitutive, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej
    REAL*8, INTENT(IN) :: qj_rel(n_vars)
    REAL*8, INTENT(IN) :: qj_org(n_vars)
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: left_matrix(n_eqns,n_vars)

    REAL*8 :: Jacob_relax(n_eqns,n_vars)
    COMPLEX*16 :: nh_terms_cmplx_impl(n_eqns)
    COMPLEX*16 :: qj_cmplx(n_vars) , qj_rel_cmplx(n_vars)

    REAL*8 :: h

    INTEGER :: i

    h = n_vars * epsilon(1.d0)
        
    ! initialize the matrix of the linearized system and the Jacobian
    
    left_matrix(1:n_eqns,1:n_vars) = 0.D0
    Jacob_relax(1:n_eqns,1:n_vars) = 0.D0
    
    ! evaluate the jacobian of the non-hyperbolic terms
    
    DO i=1,n_vars
     
       left_matrix(i,i) = coeff_f(i) * qj_org(i)
  
       IF ( implicit_flag(i) ) THEN 

          qj_rel_cmplx(1:n_vars) = qj_rel(1:n_vars)
          qj_rel_cmplx(i) = DCMPLX(qj_rel(i), h)
          
          qj_cmplx = qj_rel_cmplx * qj_org
          
          CALL eval_nonhyperbolic_terms( Bj , Bprimej , c_qj = qj_cmplx ,       &
               c_nh_term_impl = nh_terms_cmplx_impl ) 
          
          Jacob_relax(1:n_eqns,i) = coeff_f(i) *                                &
               AIMAG(nh_terms_cmplx_impl) / h
          
          left_matrix(1:n_eqns,i) = left_matrix(1:n_eqns,i) - dt * a_diag       &
               * Jacob_relax(1:n_eqns,i)
          
       END IF

    END DO

  END SUBROUTINE eval_jacobian

  !******************************************************************************
  !> \brief Evaluate the explicit terms 
  !
  !> This subroutine evaluate the explicit terms (non-fluxes) of the non-linear 
  !> system with respect to the conservative variables.
  !> \param[in]    q_expl          conservative variables 
  !> \param[out]   expl_terms      explicit terms
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_explicit_terms( q_expl , expl_terms )

    USE constitutive, ONLY : eval_explicit_forces

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells)
    REAL*8, INTENT(OUT) :: expl_terms(n_eqns,comp_cells)

    REAL*8 :: qc(n_vars)      !< conservative variables 
    REAL*8 :: expl_forces_term(n_eqns)      !< conservative variables 

    INTEGER :: j

    DO j = 1,comp_cells

       qc = q_expl(1:n_vars,j)

       CALL eval_explicit_forces( B_cent(j) , B_prime(j) , qc , expl_forces_term )

       expl_terms(1:n_eqns,j) =  expl_forces_term

    END DO
    
  END SUBROUTINE eval_explicit_terms

  !******************************************************************************
  !> \brief Semidiscrete finite volume central scheme
  !
  !> This subroutine solve the hyperbolic part of the system of the eqns,
  !> with a modified finite volume scheme from Kurganov et al. 2001, 
  !> coupled with a MUSCL-Hancock scheme (Van Leer, 1984) applied to a
  !> set of physical variables derived from the conservative vriables.
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_hyperbolic_terms( q_expl , F_x )

    ! External variables
    USE geometry, ONLY : dx
    USE parameters, ONLY : solver_scheme

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells)
    REAL*8, INTENT(OUT) :: F_x(n_eqns,comp_cells)

    REAL*8 :: q_old(n_vars,comp_cells)

    INTEGER :: i, j      !< loop counters
    
    REAL*8 :: h_new

    q_old = q

    q = q_expl

    IF ( reconstr_variables .EQ. 0 ) THEN
       
       ! Linear reconstruction of the conservative variables at the interfaces
       CALL reconstruction_cons
       
    ELSEIF ( reconstr_variables .EQ. 0 ) THEN
       
       ! Linear reconstruction of the physical var. (h+B,u) at the interfaces
       CALL reconstruction_phys1
       
    ELSE 
       
       ! Linear reconstruction of the physical var. (h,u) at the interfaces
       CALL reconstruction_phys2
       
    END IF
    
    ! Evaluation of the maximum local speeds at the interfaces
    CALL eval_speeds
    
    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )
       
    CASE ("LxF")
       
       CALL eval_flux_LxF
       
    CASE ("GFORCE")
       
       CALL eval_flux_GFORCE
       
    CASE ("KT")
       
       CALL eval_flux_KT
       
    END SELECT
    

    ! Advance in time the solution
    DO j = 1,comp_cells
          
       DO i=1,n_eqns

          F_x(i,j) = ( H_interface(i,j+1) - H_interface(i,j) ) / dx

       END DO

       h_new = q_expl(1,j) - dt * F_x(1,j) - B_cent(j)

       IF ( h_new .LT. 0.D0 ) THEN

          WRITE(*,*) 'j,h',j,h_new
          WRITE(*,*) 'dt',dt
          DO i=-2,2
             
             WRITE(*,*)  B_cent(j+i) , B_prime(j+i) , q_expl(1,j+i) , q_expl(2,j+i) 
             
          END DO

          WRITE(*,*) 'w_interface(j) ',q_interfaceL(1,j) , q_interfaceR(1,j)
          WRITE(*,*) 'w_interface(j+1) ',q_interfaceL(1,j+1) , q_interfaceR(1,j+1)

          WRITE(*,*) 'uh_interface(j) ',q_interfaceL(2,j) , q_interfaceR(2,j)
          WRITE(*,*) 'uh_interface(j+1) ',q_interfaceL(2,j+1) , q_interfaceR(2,j+1)

          WRITE(*,*) 'a(j) ',a_interfaceL(1,j) , a_interfaceR(1,j)
          WRITE(*,*) 'a(j+1) ',a_interfaceL(1,j+1) , a_interfaceR(1,j+1)

          WRITE(*,*) 'H_interface ',H_interface(i,j) , H_interface(i,j+1)
          WRITE(*,*) 
          READ(*,*) 

       END IF

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'F_x',j,F_x(:,j), H_interface(1,j+1) , H_interface(1,j)
          READ(*,*)
          
       END IF

    END DO
    
    !WRITE(*,*) 'flux left', H_interface(1,1)
    !WRITE(*,*) 'flux right', H_interface(1,comp_cells+1)
    !WRITE(*,*) 'sum fluxes',SUM(F_x(1,:))
    !READ(*,*)


    q = q_old
    
  END SUBROUTINE eval_hyperbolic_terms


  !******************************************************************************
  !> \brief Semidiscrete numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces according to Kurganov et al. 2001. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_flux_KT

    ! External procedures
    USE constitutive, ONLY : eval_fluxes

    ! External variables
    USE geometry, ONLY : dx

    IMPLICIT NONE

    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns

    REAL*8 :: flux_avg(n_eqns)   

    REAL*8 :: nc_fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: nc_fluxR(n_eqns)           !< Numerical fluxes from the eqns

    REAL*8 :: nc_flux_avg(n_eqns)

    REAL*8 :: nc_termL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: nc_termR(n_eqns)           !< Numerical fluxes from the eqns

    REAL*8 :: nc_term_avg(n_eqns)

    INTEGER :: j                      !< Loop counter
    INTEGER :: i                      !< Loop counter

    REAL*8 :: dt_loc , vel_loc

    DO j = 1 , comp_interfaces

       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceR(1:n_vars,j) ,          &
            r_flux = fluxR )

       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceL(1:n_vars,j) ,          &
            r_flux = fluxL )

       CALL average_KT( a_interfaceL(:,j) , a_interfaceR(:,j) , fluxL , fluxR , &
            flux_avg )
       
       DO i=1,n_eqns

          IF ( a_interfaceL(i,j) .EQ. a_interfaceR(i,j) ) THEN

             H_interface(i,j) = 0.D0

          ELSE

             H_interface(i,j) = flux_avg(i)                                     &
                  + ( a_interfaceR(i,j) * a_interfaceL(i,j) )                   &
                  / ( a_interfaceR(i,j) - a_interfaceL(i,j) )                   &
                  * ( q_interfaceR(i,j) - q_interfaceL(i,j) )             
             
          END IF

       END DO

!       IF ( j .EQ. comp_interfaces ) THEN
!
!          WRITE(*,*) 'flux',  H_interface(1,j) , flux_avg(1),                   &
!               ( a_interfaceR(1,j) * a_interfaceL(1,j) )                   &
!                  / ( a_interfaceR(1,j) - a_interfaceL(1,j) )                   &
!                  * ( q_interfaceR(1,j) - q_interfaceL(1,j) )
!
!       END IF

    END DO

  END SUBROUTINE eval_flux_KT

  SUBROUTINE average_KT( aL , aR , wL , wR , w_avg )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: aL(:) , aR(:)
    REAL*8, INTENT(IN) :: wL(:) , wR(:)
    REAL*8, INTENT(OUT) :: w_avg(:)

    INTEGER :: n
    INTEGER :: i 

    n = SIZE( aL )

    DO i=1,n

       IF ( aL(i) .EQ. aR(i) ) THEN

          ! w_avg(i) = 0.5D0 * ( wL(i) + wR(i) )
          w_avg(i) = 0.0D0

       ELSE

          w_avg(i) = ( aR(i) * wL(i) - aL(i) * wR(i) ) / ( aR(i) - aL(i) )  

       END IF

    END DO

  END SUBROUTINE average_KT
  
  !******************************************************************************
  !> \brief Numerical fluxes GFORCE
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************
  
  SUBROUTINE eval_flux_GFORCE
    
    ! External procedures
    USE constitutive, ONLY : eval_fluxes
    
    ! External variables
    USE geometry, ONLY : dx
    
    IMPLICIT NONE
    
    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL*8 :: flux_lf(n_eqns)         
    REAL*8 :: flux_lw(n_eqns)         
    REAL*8 :: q_lw(n_vars)           
    
    INTEGER :: j                      !< Loop counter on the interfaces
    INTEGER :: i                      !< Loop counter on the equations

    REAL*8 :: cfl_loc
    REAL*8 :: vel_loc
    REAL*8 :: dt_loc
    REAL*8 :: omega 
    
    cfl_loc = 0.9

    omega = 1.D0 / ( 1.D0 + cfl_loc )

    DO j = 1,comp_interfaces
       
       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceL(1:n_vars,j) ,          &
            r_flux = fluxL )

       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceR(1:n_vars,j) ,          &
            r_flux = fluxR )
       
       vel_loc = MAX ( ABS( a_interfaceL(1,j) ) , ABS( a_interfaceR(1,j) ) )
       
       dt_loc = cfl_loc * dx / vel_loc
       
       DO i=1,n_vars
          
          flux_lf(i) = 0.5D0 * ( fluxL(i) + fluxR(i) )  - 0.5D0 * dx /          &
               dt_loc * ( q_interfaceR(i,j) - q_interfaceL(i,j) ) 
          
          q_lw(i) = 0.5D0 * ( q_interfaceR(i,j) + q_interfaceL(i,j) )           &
               - 0.5D0 * dt_loc / dx * ( fluxR(i) - fluxL(i) )
      
       END DO
       
       CALL eval_fluxes ( B_stag(j) , r_qj = q_lw, r_flux = flux_lw )
       
       DO i=1,n_eqns

          H_interface(i,j) = ( 1.D0 - omega ) * flux_lf(i) + omega * flux_lw(i)
          
       END DO
       
    END DO

  END SUBROUTINE eval_flux_GFORCE

  !*****************************************************************************
  !> \brief Numerical fluxes Lax-Friedrichs
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !*****************************************************************************
  
  SUBROUTINE eval_flux_LxF
    
    ! External procedures
    USE constitutive, ONLY : eval_fluxes
    
    ! External variables
    USE geometry, ONLY : dx
    
    IMPLICIT NONE

    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL*8 :: flux_lf(n_eqns)         

    INTEGER :: j                      !< Loop counter
    INTEGER :: i                      !< Loop counter
    REAL*8 :: vel_loc                 !< Largest local speed
    REAL*8 :: dt_loc                  !

    DO j = 1,comp_interfaces

       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceL(1:n_eqns,j) ,          &
            r_flux = fluxL )

       CALL eval_fluxes( B_stag(j) , r_qj = q_interfaceR(1:n_eqns,j) ,          &
            r_flux = fluxR )
       
       vel_loc = MAX( ABS( a_interfaceL(1,j) ) , ABS( a_interfaceR(1,j) ) )

       dt_loc = 0.9D0 * dx / vel_loc
       
       DO i = 1,n_eqns
          
          flux_lf(i) = 0.5D0 * ( fluxR(i) + fluxL(i) ) - 0.5D0 * dx /         &
               dt_loc *( q_interfaceR(i,j) - q_interfaceL(i,j) )
          
       END DO

       DO i=1,n_eqns

          H_interface(i,j) = flux_lf(i)
          
       END DO
       
    END DO
    
  END SUBROUTINE eval_flux_LxF

  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of conservative variables describing the state of the
  !> system. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 17/11/2016
  !******************************************************************************

  SUBROUTINE reconstruction_cons

    ! External procedures
    USE constitutive, ONLY : qc_to_qp , qp_to_qc
    USE parameters, ONLY : limiter

    ! External variables
    USE geometry, ONLY : x_comp , x_stag
    USE parameters, ONLY : reconstr_coeff

    IMPLICIT NONE

    REAL*8 :: qc(n_vars)      !< conservative variables
    REAL*8 :: qcL(n_vars)     !< conservative variables at the left edge of the cells
    REAL*8 :: qcR(n_vars)     !< conservative variables at the rightedge of the cells

    REAL*8 :: qc_stencil(3)   !< conservative variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: qc_prime        !< conservative variables slope

    REAL*8 :: qp_check(n_vars)

    REAL*8 :: dxL , dxR

    INTEGER :: j              !< loop counter
    INTEGER :: i              !< loop counter

    DO i=1,n_vars

       ! Slope in the first cell
       IF ( bcL(i)%flag .EQ. 0 ) THEN
          ! Dirichelet boundary condition

          qc_stencil(1) = bcL(i)%value
          qc_stencil(2:3) = q(i,1:2)
          
          x_stencil(1) = x_stag(1)
          x_stencil(2:3) = x_comp(1:2)
          
          CALL limit( qc_stencil , x_stencil , limiter(i) , qc_prime ) 

       ELSEIF ( bcL(i)%flag .EQ. 1 ) THEN
          ! Neumann boundary condition (fixed slope)
        
          qc_prime = bcL(i)%value

       ELSEIF ( bcL(i)%flag .EQ. 2 ) THEN
          ! Linear extrapolation from internal values

          qc_prime = ( q(i,2) - q(i,1) ) / ( x_comp(2) - x_comp(1) )
          
       END IF

       ! Linear reconstruction of the conservative variables at the boundaries
       ! of for the first cell
       dxL = x_comp(1) - x_stag(1)
       qcL(i) = q(i,1) - reconstr_coeff * dxL * qc_prime

       dxR = x_stag(2) - x_comp(1)
       qcR(i) = q(i,1) + reconstr_coeff * dxR * qc_prime

       ! Positivity preserving reconstruction for h
       IF (i.eq.1) THEN
          
          IF ( qcR(i) .LT. B_stag(2) ) THEN
             
             qc_prime = ( B_stag(2) - q(i,1) ) / dxR
             
             qcL(i) = q(i,1) - reconstr_coeff * dxL * qc_prime
             
             qcR(i) = q(i,1) + reconstr_coeff * dxR * qc_prime
             
          ENDIF
          
          IF ( qcL(i) .LT. B_stag(1) ) THEN
             
             qc_prime = ( q(i,1) - B_stag(1) ) / dxL
             
             qcL(i) = q(i,1) - reconstr_coeff * dxL * qc_prime
             
             qcR(i) = q(i,1) + reconstr_coeff * dxR * qc_prime
             
          ENDIF
          
       ENDIF
       
    END DO

    ! Correction for small values (Eqs. 2.17 and 2.21 K&P)
    CALL qc_to_qp( qcL , B_stag(1) , qp_check )
    CALL qp_to_qc( qp_check , B_stag(1) , qcL )

    CALL qc_to_qp( qcR , B_stag(2) , qp_check )
    CALL qp_to_qc( qp_check , B_stag(2) , qcR )
 
    
    ! Value at the right of the first interface
    q_interfaceR(1:n_vars,1) = qcL
    ! Value at the left of the second interface
    q_interfaceL(1:n_vars,2) = qcR

    ! Values at the left of the first interface (out of the domain)
    DO i=1,n_vars

       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          q_interfaceR(i,1) = bcL(i)%value 
          q_interfaceL(i,1) = bcL(i)%value 
        
       ELSEIF ( bcL(i)%flag .EQ. 1 ) THEN

          q_interfaceL(i,1) = q_interfaceR(i,1) 
  
       ELSE
          
          q_interfaceL(i,1) = qcL(i)
          
       END IF

    END DO

    ! Linear reconstruction of the physical variables for the internal cells
    DO j = 2,comp_cells-1
       
       x_stencil(1:3) = x_comp(j-1:j+1)
       
       DO i=1,n_vars
          
          qc_stencil = q(i,j-1:j+1)
          
          CALL limit( qc_stencil , x_stencil , limiter(i) , qc_prime ) 
          
          ! Linear reconstruction at the left bdry of the cell
          dxL = x_comp(j) - x_stag(j)
          qcL(i) = q(i,j) - reconstr_coeff * dxL * qc_prime
          
          ! Linear reconstruction at the right bdry of the cell
          dxR = x_stag(j+1) - x_comp(j)
          qcR(i) = q(i,j) + reconstr_coeff * dxR * qc_prime
          
          ! Positivity preserving reconstruction for h
          IF ( i .EQ. 1 ) THEN
             
             IF ( qcR(i) .LT. B_stag(j+1) ) THEN
                
                qc_prime = ( B_stag(j+1) - q(i,j) ) / dxR
                
                qcL(i) = q(i,j) - reconstr_coeff * dxL * qc_prime
                
                qcR(i) = q(i,j) + reconstr_coeff * dxR * qc_prime
                
             ENDIF
             
             IF ( qcL(i) .LT. B_stag(j) ) THEN
                
                qc_prime = ( q(i,j) - B_stag(j) ) / dxL
                
                qcL(i) = q(i,j) - reconstr_coeff * dxL * qc_prime
                
                qcR(i) = q(i,j) + reconstr_coeff * dxR * qc_prime
                
             ENDIF
             
          ENDIF
          
       END DO
       
              
       ! Correction for small values (Eqs. 2.17 and 2.21 K&P)
       CALL qc_to_qp( qcL , B_stag(j) , qp_check )
       CALL qp_to_qc( qp_check , B_stag(j) , qcL )
       
 
       CALL qc_to_qp( qcR , B_stag(j+1) , qp_check )
       CALL qp_to_qc( qp_check , B_stag(j+1) , qcR )
       
       ! Value at the right of the j-th interface
       q_interfaceR(:,j) = qcL
       ! Value at the left of the j+1-th interface
       q_interfaceL(:,j+1) = qcR

    END DO

    ! Linear reconstruction of the physical variables at the boundaries
    ! of the last cell (j=comp_cells)

    DO i=1,n_vars
       
       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qc_stencil(1:2) = q(i,comp_cells-1:comp_cells)
          qc_stencil(3) = bcR(i)%value
          
          x_stencil(1:2) = x_comp(comp_cells-1:comp_cells)
          x_stencil(3) = x_stag(comp_interfaces)
          
          CALL limit( qc_stencil , x_stencil , limiter(i) , qc_prime ) 

       ELSEIF ( bcR(i)%flag .EQ. 1 ) THEN
        
          qc_prime = bcR(i)%value
          
       ELSEIF ( bcR(i)%flag .EQ. 2 ) THEN
        
          qc_prime = ( q(i,comp_cells) - q(i,comp_cells-1) ) /                &
               ( x_comp(comp_cells) - x_comp(comp_cells-1) )

       END IF

       dxL = x_comp(comp_cells) - x_stag(comp_interfaces-1)
       qcL(i) = q(i,comp_cells) - reconstr_coeff * dxL * qc_prime
       
       dxR = x_stag(comp_interfaces) - x_comp(comp_cells)
       qcR(i) = q(i,comp_cells) + reconstr_coeff * dxR * qc_prime

       ! Positivity preserving reconstruction for h
       IF (i.eq.1) THEN
          
          IF ( qcR(i) .LT. B_stag(comp_interfaces) ) THEN
             
             qc_prime = ( B_stag(comp_interfaces) - q(i,comp_cells) ) / dxR
             
             qcL(i) = q(i,comp_cells) - reconstr_coeff * dxL * qc_prime
             
             qcR(i) = q(i,comp_cells) + reconstr_coeff * dxR * qc_prime
             
          ENDIF
          
          IF ( qcL(i) .LT. B_stag(comp_interfaces-1) ) THEN
             
             qc_prime = ( q(i,comp_cells) - B_stag(comp_cells) ) / dxL
             
             qcL(i) = q(i,comp_cells) - reconstr_coeff * dxL * qc_prime
             
             qcR(i) = q(i,comp_cells) + reconstr_coeff * dxR * qc_prime
             
          ENDIF
          
       ENDIF
       
    END DO
    
    ! Correction for small values (Eqs. 2.17 and 2.21 K&P)
    CALL qc_to_qp( qcL , B_stag(comp_cells) , qp_check )
    CALL qp_to_qc( qp_check , B_stag(comp_cells) , qcL )
    
    CALL qc_to_qp( qcR , B_stag(comp_cells+1) , qp_check )
    CALL qp_to_qc( qp_check , B_stag(comp_cells+1) , qcR )
    
    q_interfaceR(:,comp_interfaces-1) = qcL
    q_interfaceL(:,comp_interfaces) = qcR
  
    ! Values at the right of the last interface (out of the domain)
    DO i=1,n_vars

       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          q_interfaceL(i,comp_interfaces) = bcR(i)%value 
          q_interfaceR(i,comp_interfaces) = bcR(i)%value 
          
       ELSEIF ( bcR(i)%flag .EQ. 1 ) THEN
          
          q_interfaceR(i,comp_interfaces) = q_interfaceL(i,comp_interfaces) 
          
       ELSE
          
          q_interfaceR(i,comp_interfaces) = qcR(i)
          
       END IF

    END DO

  END SUBROUTINE reconstruction_cons

  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of physical variables (h+B,u) describing the state of the
  !> system. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 17/11/2016
  !******************************************************************************


  SUBROUTINE reconstruction_phys1

    ! External procedures
    USE constitutive, ONLY : qc_to_qp , qp_to_qc
    USE constitutive, ONLY : qc_to_qp2 , qp2_to_qc
    USE parameters, ONLY : limiter

    ! External variables
    USE geometry, ONLY : x_comp , x_stag
    USE parameters, ONLY : reconstr_coeff

    IMPLICIT NONE

    REAL*8 :: qc(n_vars)      !< conservative variables
    REAL*8 :: qpL(n_vars)     !< physical variables at the left edge of the cells
    REAL*8 :: qpR(n_vars)     !< physical variables at the rightedge of the cells
    REAL*8 :: qp_bdry(n_vars) !< physical variables outside the domain

    REAL*8 :: qp_stencil(3)   !< physical variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: qp_prime        !< physical variables slope

    REAL*8 :: dxL , dxR

    INTEGER :: j              !< loop counter
    INTEGER :: i              !< loop counter


    ! Convert the conservative variables to the physical variables
    DO j = 1,comp_cells

       qc = q(1:n_vars,j)

       CALL qc_to_qp( qc , B_cent(j) , qp(1:n_vars,j) )

    END DO

    DO i=1,n_vars

       ! Slope in the first cell
       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qp_stencil(1) = bcL(i)%value
          qp_stencil(2:3) = qp(i,1:2)
          
          x_stencil(1) = x_stag(1)
          x_stencil(2:3) = x_comp(1:2)
          
          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

       ELSEIF ( bcL(i)%flag .EQ. 1 ) THEN
        
          qp_prime = bcL(i)%value

       ELSEIF ( bcL(i)%flag .EQ. 2 ) THEN
        
          qp_prime = ( qp(i,2) - qp(i,1) ) / ( x_comp(2) - x_comp(1) )

       END IF

       ! Linear reconstruction of the physical variables at the boundaries
       ! of for the first cell
       dxL = x_comp(1) - x_stag(1)
       qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

       dxR = x_stag(2) - x_comp(1)
       qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

       ! positivity preserving reconstruction for h
       IF(i.eq.1)THEN

         IF(qpR(i).LT.B_stag(2))THEN

            qp_prime=(B_stag(2)-qp(i,1))/dxR

            qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

            qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

         ENDIF

         IF(qpL(i).LT.B_stag(1))THEN

            qp_prime=(qp(i,1)-B_stag(1))/dxL

            qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

            qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

         ENDIF

       ENDIF
       
    END DO
    
    DO i=1,n_vars
       
       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qpL(i) = bcL(i)%value
          
       END IF
       
    END DO


    ! Convert back from physical to conservative variables
    CALL qp_to_qc( qpL , B_stag(1) , q_interfaceR(1:n_vars,1) )
    CALL qp_to_qc( qpR , B_stag(2) , q_interfaceL(1:n_vars,2) )


    DO i=1,n_vars

       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qp_bdry(i) = bcL(i)%value 
          
       ELSE
          
          !qp_bdry(i) = qpL(i)

          ! fixed wall
          IF(i.eq.2)THEN

             qp_bdry(i) = -qpL(i)

          ELSE

             qp_bdry(i) = qpL(i)

          ENDIF
          
       END IF

    END DO

    CALL qp_to_qc( qp_bdry , B_stag(1) , q_interfaceL(:,1) )

    ! Linear reconstruction of the physical variables for the internal cells
    DO j = 2,comp_cells-1

       x_stencil(1:3) = x_comp(j-1:j+1)

       DO i=1,n_vars

          qp_stencil = qp(i,j-1:j+1)

          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

          dxL = x_comp(j) - x_stag(j)
          qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime
          
          dxR = x_stag(j+1) - x_comp(j)
          qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

          ! positivity preserving reconstruction for h
          IF(i.eq.1)THEN
 
            IF(qpR(i).LT.B_stag(j+1))THEN
 
              qp_prime=(B_stag(j+1)-qp(i,j))/dxR

              qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime

              qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

            ENDIF

            IF(qpL(i).LT.B_stag(j))THEN

              qp_prime=(qp(i,j)-B_stag(j))/dxL

              qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime

              qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

            ENDIF

          ENDIF
          
       END DO

       ! Convert back from physical to conservative variables
       CALL qp_to_qc( qpL , B_stag(j) , q_interfaceR(:,j) )
       CALL qp_to_qc( qpR , B_stag(j+1) , q_interfaceL(:,j+1) )

    END DO

    ! Linear reconstruction of the physical variables at the boundaries
    ! of the last cell (j=comp_cells)

    DO i=1,n_vars
       
       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qp_stencil(1:2) = qp(i,comp_cells-1:comp_cells)
          qp_stencil(3) = bcR(i)%value
          
          x_stencil(1:2) = x_comp(comp_cells-1:comp_cells)
          x_stencil(3) = x_stag(comp_interfaces)
          
          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

       ELSEIF ( bcR(i)%flag .EQ. 1 ) THEN
        
          qp_prime = bcR(i)%value
          
       ELSEIF ( bcR(i)%flag .EQ. 2 ) THEN
        
          qp_prime = ( qp(i,comp_cells) - qp(i,comp_cells-1) ) /                &
               ( x_comp(comp_cells) - x_comp(comp_cells-1) )

       END IF

       dxL = x_comp(comp_cells) - x_stag(comp_interfaces-1)
       qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime
       
       dxR = x_stag(comp_interfaces) - x_comp(comp_cells)
       qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

       ! positivity preserving reconstruction for h
       IF(i.eq.1)THEN

         IF(qpR(i).LT.B_stag(comp_interfaces))THEN

           qp_prime=(B_stag(comp_interfaces)-qp(i,comp_cells))/dxR

           qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime

           qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

         ENDIF

         IF(qpL(i).LT.B_stag(comp_interfaces-1))THEN

           qp_prime=(qp(i,comp_cells)-B_stag(comp_cells))/dxL

           qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime

           qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

         ENDIF

       ENDIF
       
    END DO
    
    DO i=1,n_vars

       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qpR(i) = bcR(i)%value
          
       END IF

    END DO

    ! Convert back from physical to conservative variables
    CALL qp_to_qc( qpL , B_stag(comp_interfaces-1) , q_interfaceR(:,comp_interfaces-1) )
    CALL qp_to_qc( qpR , B_stag(comp_interfaces) , q_interfaceL(:,comp_interfaces) )
  
    DO i=1,n_vars

       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qp_bdry(i) = bcR(i)%value 
          
       ELSE
          
          !qp_bdry(i) = qpR(i)

          ! fixed wall
          IF(i.eq.2)THEN

             qp_bdry(i) = -qpR(i)

          ELSE

             qp_bdry(i) = qpR(i)

          ENDIF
          
       END IF

    END DO

    CALL qp_to_qc( qp_bdry , B_stag(comp_interfaces) ,                          &
        q_interfaceR(:,comp_interfaces) )

  END SUBROUTINE reconstruction_phys1

  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of physical variables (h,u) describing the state of the
  !> system. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 17/11/2016
  !******************************************************************************

  SUBROUTINE reconstruction_phys2

    ! External procedures
    USE constitutive, ONLY : qc_to_qp , qp_to_qc
    USE constitutive, ONLY : qc_to_qp2 , qp2_to_qc
    USE parameters, ONLY : limiter

    ! External variables
    USE geometry, ONLY : x_comp , x_stag
    USE parameters, ONLY : reconstr_coeff

    IMPLICIT NONE

    REAL*8 :: qc(n_vars)      !< conservative variables
    REAL*8 :: qpL(n_vars)     !< physical variables at the left edge of the cells
    REAL*8 :: qpR(n_vars)     !< physical variables at the rightedge of the cells
    REAL*8 :: qp_bdry(n_vars) !< physical variables outside the domain

    REAL*8 :: qp_stencil(3)   !< physical variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: qp_prime        !< physical variables slope

    REAL*8 :: dxL , dxR

    INTEGER :: j              !< loop counter
    INTEGER :: i              !< loop counter


    ! Convert the conservative variables to the physical variables
    DO j = 1,comp_cells

       qc = q(1:n_vars,j)

       CALL qc_to_qp2( qc , B_cent(j) , qp(1:n_vars,j) )

    END DO

    DO i=1,n_vars

       ! Slope in the first cell
       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qp_stencil(1) = bcL(i)%value
          qp_stencil(2:3) = qp(i,1:2)
          
          x_stencil(1) = x_stag(1)
          x_stencil(2:3) = x_comp(1:2)
          
          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

       ELSEIF ( bcL(i)%flag .EQ. 1 ) THEN
        
          qp_prime = bcL(i)%value

       ELSEIF ( bcL(i)%flag .EQ. 2 ) THEN
        
          qp_prime = ( qp(i,2) - qp(i,1) ) / ( x_comp(2) - x_comp(1) )

       END IF

       ! Linear reconstruction of the physical variables at the boundaries
       ! of for the first cell
       dxL = x_comp(1) - x_stag(1)
       qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

       dxR = x_stag(2) - x_comp(1)
       qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

       ! positivity preserving reconstruction for h
       IF ( i .EQ. 1 ) THEN

         IF ( qpR(i) .LT. 0.D0 ) THEN

            qp_prime = - qp(i,1) / dxR

            qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

            qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

         ENDIF

         IF ( qpL(i) .LT. 0.D0 ) THEN

            qp_prime = qp(i,1) / dxL

            qpL(i) = qp(i,1) - reconstr_coeff * dxL * qp_prime

            qpR(i) = qp(i,1) + reconstr_coeff * dxR * qp_prime

         ENDIF

       ENDIF
       
    END DO
    
    DO i=1,n_vars
       
       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qpL(i) = bcL(i)%value
          
       END IF
       
    END DO


    ! Convert back from physical to conservative variables
    CALL qp2_to_qc( qpL , B_stag(1) , q_interfaceR(1:n_vars,1) )
    CALL qp2_to_qc( qpR , B_stag(2) , q_interfaceL(1:n_vars,2) )


    DO i=1,n_vars

       IF ( bcL(i)%flag .EQ. 0 ) THEN
          
          qp_bdry(i) = bcL(i)%value 
          
       ELSE
          
          !qp_bdry(i) = qpL(i)

          ! fixed wall
          IF(i.eq.2)THEN

             qp_bdry(i) = -qpL(i)

          ELSE

             qp_bdry(i) = qpL(i)

          ENDIF
          
       END IF

    END DO

    CALL qp2_to_qc( qp_bdry , B_stag(1) , q_interfaceL(:,1) )

    ! Linear reconstruction of the physical variables for the internal cells
    DO j = 2,comp_cells-1

       x_stencil(1:3) = x_comp(j-1:j+1)

       DO i=1,n_vars

          qp_stencil = qp(i,j-1:j+1)

          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

          dxL = x_comp(j) - x_stag(j)
          qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime
          
          dxR = x_stag(j+1) - x_comp(j)
          qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

          ! positivity preserving reconstruction for h
          IF ( i .EQ. 1 ) THEN
 
            IF ( qpR(i) .LT. 0.D0 ) THEN
 
              qp_prime = - qp(i,j) / dxR

              qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime

              qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

            ENDIF

            IF ( qpL(i) .LT. 0.D0 ) THEN

              qp_prime = qp(i,j) / dxL

              qpL(i) = qp(i,j) - reconstr_coeff * dxL * qp_prime

              qpR(i) = qp(i,j) + reconstr_coeff * dxR * qp_prime

            ENDIF

          ENDIF
          
       END DO

       ! Convert back from physical to conservative variables
       CALL qp2_to_qc( qpL , B_stag(j) , q_interfaceR(:,j) )
       CALL qp2_to_qc( qpR , B_stag(j+1) , q_interfaceL(:,j+1) )

    END DO

    ! Linear reconstruction of the physical variables at the boundaries
    ! of the last cell (j=comp_cells)

    DO i=1,n_vars
       
       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qp_stencil(1:2) = qp(i,comp_cells-1:comp_cells)
          qp_stencil(3) = bcR(i)%value
          
          x_stencil(1:2) = x_comp(comp_cells-1:comp_cells)
          x_stencil(3) = x_stag(comp_interfaces)
          
          CALL limit( qp_stencil , x_stencil , limiter(i) , qp_prime ) 

       ELSEIF ( bcR(i)%flag .EQ. 1 ) THEN
        
          qp_prime = bcR(i)%value
          
       ELSEIF ( bcR(i)%flag .EQ. 2 ) THEN
        
          qp_prime = ( qp(i,comp_cells) - qp(i,comp_cells-1) ) /                &
               ( x_comp(comp_cells) - x_comp(comp_cells-1) )

       END IF

       dxL = x_comp(comp_cells) - x_stag(comp_interfaces-1)
       qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime
       
       dxR = x_stag(comp_interfaces) - x_comp(comp_cells)
       qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

       ! positivity preserving reconstruction for h
       IF ( i .EQ. 1 ) THEN

         IF ( qpR(i) .LT. 0.D0 ) THEN

           qp_prime = - qp(i,comp_cells) / dxR

           qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime

           qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

         ENDIF

         IF ( qpL(i) .LT. 0.D0 ) THEN

           qp_prime = qp(i,comp_cells) / dxL

           qpL(i) = qp(i,comp_cells) - reconstr_coeff * dxL * qp_prime

           qpR(i) = qp(i,comp_cells) + reconstr_coeff * dxR * qp_prime

         ENDIF

       ENDIF
       
    END DO
    
    DO i=1,n_vars

       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qpR(i) = bcR(i)%value
          
       END IF

    END DO

    ! Convert back from physical to conservative variables
    CALL qp2_to_qc( qpL , B_stag(comp_interfaces-1) , q_interfaceR(:,comp_interfaces-1) )
    CALL qp2_to_qc( qpR , B_stag(comp_interfaces) , q_interfaceL(:,comp_interfaces) )
  
    DO i=1,n_vars

       IF ( bcR(i)%flag .EQ. 0 ) THEN
          
          qp_bdry(i) = bcR(i)%value 
          
       ELSE
          
          !qp_bdry(i) = qpR(i)

          ! fixed wall
          IF(i.eq.2)THEN

             qp_bdry(i) = -qpR(i)

          ELSE

             qp_bdry(i) = qpR(i)

          ENDIF
          
       END IF

    END DO

    CALL qp2_to_qc( qp_bdry , B_stag(comp_interfaces) ,                          &
        q_interfaceR(:,comp_interfaces) )

  END SUBROUTINE reconstruction_phys2

  !******************************************************************************
  !> \brief Characteristic speeds
  !
  !> This subroutine evaluates the largest characteristic speed at the
  !> cells interfaces from the reconstructed states.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_speeds

    ! External procedures
    USE constitutive, ONLY : eval_local_speeds2

    IMPLICIT NONE

    REAL*8 :: abslambdaL_min(n_vars) , abslambdaL_max(n_vars)
    REAL*8 :: abslambdaR_min(n_vars) , abslambdaR_max(n_vars)
    REAL*8 :: min_r(n_vars) , max_r(n_vars)

    INTEGER :: j

    DO j = 1 , comp_interfaces

       CALL eval_local_speeds2( q_interfaceL(:,j) , B_stag(j) , abslambdaL_min , &
               abslambdaL_max )
       CALL eval_local_speeds2( q_interfaceR(:,j) , B_stag(j) , abslambdaR_min , &
               abslambdaR_max )

       min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0D0)
       max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0D0)
       
       a_interfaceL(:,j) = min_r
       a_interfaceR(:,j) = max_r

    END DO

  END SUBROUTINE eval_speeds



  !******************************************************************************
  !> \brief Slope limiter
  !
  !> This subroutine limits the slope of the linear reconstruction of 
  !> the physical variables, accordingly to the parameter "solve_limiter":\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  !> \param[in]     v             3-point stencil value array 
  !> \param[in]     z             3-point stencil location array 
  !> \param[out]    slope_lim     limited slope         
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE limit( v , z , limiter , slope_lim )

    USE parameters, ONLY : theta

    IMPLICIT none

    REAL*8, INTENT(IN) :: v(3)
    REAL*8, INTENT(IN) :: z(3)
    INTEGER, INTENT(IN) :: limiter

    REAL*8, INTENT(OUT) :: slope_lim

    REAL*8 :: a , b , c

    REAL*8 :: sigma1 , sigma2
 
    a = ( v(3) - v(2) ) / ( z(3) - z(2) )
    b = ( v(2) - v(1) ) / ( z(2) - z(1) )
    c = ( v(3) - v(1) ) / ( z(3) - z(1) )

    SELECT CASE (limiter)

    CASE ( 0 )

       slope_lim = 0.D0

    CASE ( 1 )

       ! minmod

       slope_lim = minmod(a,b)

    CASE ( 2 )

       ! superbee

       sigma1 = minmod( a , 2.D0*b )
       sigma2 = minmod( 2.D0*a , b )
       slope_lim = maxmod( sigma1 , sigma2 )

    CASE ( 3 )

       ! van_leer

       slope_lim = minmod( c , theta * minmod( a , b ) )

    END SELECT

  END SUBROUTINE limit


  REAL*8 FUNCTION minmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .LE. 0.D0 ) THEN

       minmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       minmod = 0.5 * ( sa+sb ) * MIN( ABS(a) , ABS(b) )

    END IF

  END FUNCTION minmod

  REAL*8 function maxmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .EQ. 0.d0 ) THEN

       maxmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       maxmod = 0.5 * ( sa+sb ) * MAX( ABS(a) , ABS(b) )

    END IF

  END function maxmod

END MODULE solver
