!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive

  USE parameters, ONLY : n_eqns , n_vars

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name

  !--------- Constants for the equations of state -----------------------------


  COMPLEX*16 :: h      !< height
  COMPLEX*16 :: u      !< velocity (x direction)

  !> gravitational acceleration
  REAL*8 :: grav

CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters, ONLY : n_nh

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.

    n_nh = COUNT( implicit_flag )

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$\alpha_i, u_i, \rho_i, x_i, \f$)
  !> for the two phases and the mixture denisity \f$ \rho \f$.
  !> \param[in]    r_qj     real conservative variables 
  !> \param[in]    c_qj     complex conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE phys_var(Bj,r_qj,c_qj)
    
    USE COMPLEXIFY
    USE parameters, ONLY : eps_sing
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16 :: qj(n_vars)

    IF ( present(c_qj) ) THEN

       qj = c_qj

    ELSE

       qj = DCMPLX(r_qj)

    END IF

    h = qj(1) - DCMPLX( Bj , 0.D0 )

    
    ! Correction for small values of h (Eq. 2.17 K&P 2007)
    IF ( REAL( h ) .GT. eps_sing ** 0.25D0 ) THEN

       u = qj(2) / h 

    ELSE

       u = DSQRT(2.D0) * h * qj(2) / CDSQRT( h**4 + DCMPLX(eps_sing,0.D0) )

    END IF

  END SUBROUTINE phys_var

  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on u and h, changing
  !> u when h is small (Eq. 2.17 K&P 2007).
  !> \date 10/04/2012
  !******************************************************************************

  SUBROUTINE eval_local_speeds(qj,Bj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , u_temp

    CALL phys_var(Bj,r_qj = qj)
    
    vel_min(1:n_eqns) = REAL(u) - DSQRT( grav * REAL(h) ) 
    vel_max(1:n_eqns) = REAL(u) + DSQRT( grav * REAL(h) )

  END SUBROUTINE eval_local_speeds
  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on u and h.
  !> \date 10/04/2012
  !******************************************************************************

  SUBROUTINE eval_local_speeds2(qj,Bj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , u_temp

    h_temp = qj(1) - Bj

    IF ( h_temp .NE. 0.D0 ) THEN

       u_temp = qj(2) / h_temp

    ELSE

       u_temp = 0.D0

    END IF

    vel_min(1:n_eqns) = REAL(u_temp) - DSQRT( grav * REAL(h_temp) ) 
    vel_max(1:n_eqns) = REAL(u_temp) + DSQRT( grav * REAL(h_temp) )


  END SUBROUTINE eval_local_speeds2



  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h+B \f$
  !> - qp(2) = \f$ u \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces (limiters are applied to the reconstructed slopes) and for 
  !> the boundary conditions. 
  !> Please note that the physical variable u is modified according to Eq. 2.17 
  !> of K&P 2007 by the call to the subroutine phys_var.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(B,r_qj = qc)
          
    qp(1) = REAL(h+B)
    qp(2) = REAL(u)

    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h + B \f$
  !> - qp(2) = \f$ hu \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_hB      !> batimetry + height 
    REAL*8 :: r_u       !> velocity
    
    r_hB = qp(1)
    r_u = qp(2)
    
    qc(1) = r_hB
    qc(2) = ( r_hB - B ) * r_u

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Conservative to 2nd set of physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ u \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp2(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(B,r_qj = qc)
          
    qp(1) = REAL(h)
    qp(2) = REAL(u)

    
  END SUBROUTINE qc_to_qp2


  !******************************************************************************
  !> \brief From 2nd set of physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ u \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp2_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_h      !> batimetry + height 
    REAL*8 :: r_u       !> velocity
    
    r_h = qp(1)
    r_u = qp(2)
    
    qc(1) = r_h + B
    qc(2) = r_h * r_u

  END SUBROUTINE qp2_to_qc

  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qj, accordingly to the equations for the single temperature
  !> model introduced in Romenki et al. 2010.
  !> \date 22/11/2016
  !> \param[in]     c_qj     complex conservative variables 
  !> \param[in]     r_qj     real conservative variables 
  !> \param[out]    c_flux   complex analytical fluxes    
  !> \param[out]    r_flux   real analytical fluxes    
  !******************************************************************************
  
  SUBROUTINE eval_fluxes(Bj,c_qj,r_qj,c_flux,r_flux)
    
    USE COMPLEXify 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: Bj
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_flux(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_flux(n_eqns)
    
    COMPLEX*16 :: qj(n_vars)
    COMPLEX*16 :: flux(n_eqns)

    INTEGER :: i 
    
    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF
   
    flux(1) = qj(2)
    
    ! check if h is small
    IF((qj(1)-Bj).GT.10.d0**(-8))THEN

      flux(2) = qj(2)**2.D0 / ( qj(1) - Bj ) + 0.5D0 * grav * ( qj(1) - Bj ) ** 2.D0 

    ELSE

      flux(2) = 0.d0

    ENDIF    

    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       c_flux = flux

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       r_flux = REAL( flux )

    END IF
  
  END SUBROUTINE eval_fluxes

  !******************************************************************************
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms( Bj , Bprimej , c_qj , c_nh_term_impl ,   &
       r_qj , r_nh_term_impl )

    USE COMPLEXIFY 
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: nh_term(n_eqns)

    COMPLEX*16 :: relaxation_term(n_eqns)
    COMPLEX*16 :: heat , drag
    COMPLEX*16 :: vel_term

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i
    
    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the relaxation terms
    relaxation_term(1:n_eqns) = DCMPLX(0.D0) 

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0)

    nh_term = relaxation_term + forces_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF 

  END SUBROUTINE eval_nonhyperbolic_terms

  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the forces to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity)
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    expl_forces_term   forces term
  !******************************************************************************

  SUBROUTINE eval_explicit_forces( Bj , Bprimej , qj , expl_forces_term )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej

    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    REAL*8, INTENT(OUT) :: expl_forces_term(n_eqns)  !< explicit forces 

    expl_forces_term(1:n_eqns) = 0.D0

    CALL phys_var(Bj,r_qj = qj)

    expl_forces_term(2) = grav * h * Bprimej

    ! friction term
    expl_forces_term(2) = expl_forces_term(2) + ( 0.001 / (1.D0+10.D0*h) ) * u

  END SUBROUTINE eval_explicit_forces

END MODULE constitutive

    
