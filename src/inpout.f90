!********************************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************

MODULE inpout

  ! -- Variables for the namelist RUN_PARAMETERS
  USE parameters, ONLY : solver_scheme, max_dt , t_start , t_end ,              &
       dt_output , cfl, limiter , theta, reconstr_coeff ,                       &
       interfaces_relaxation , n_RK

  USE solver, ONLY : verbose_level

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry, ONLY : x0 , xN , comp_cells
  USE geometry, ONLY : batimetry_profile , n_batimetry_profile
  USE init, ONLY : riemann_interface

  ! -- Variables for the namelist LEFT_STATE
  USE init, ONLY : hB_L , u_L

  ! -- Variables for the namelist RIGHT_STATE
  USE init, ONLY : hB_R , u_R

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONS
  USE parameters, ONLY : bc

  ! -- Variables for the namelist SOURCE_PARAMETERS
  USE constitutive, ONLY : grav
  USE constitutive, ONLY : friction_model

  IMPLICIT NONE

  CHARACTER(LEN=40) :: run_name           !< Name of the run
  CHARACTER(LEN=40) :: bak_name           !< Backup file for the parameters
  CHARACTER(LEN=40) :: input_file         !< File with the run parameters
  CHARACTER(LEN=40) :: output_file        !< Name of the output files
  CHARACTER(LEN=40) :: restart_file       !< Name of the restart file 
  CHARACTER(LEN=40) :: dakota_file        !< Name of the dakota file 

  INTEGER, PARAMETER :: input_unit = 7    !< Input data unit
  INTEGER, PARAMETER :: backup_unit = 8   !< Backup input data unit
  INTEGER, PARAMETER :: output_unit = 9   !< Output data unit
  INTEGER, PARAMETER :: restart_unit = 10 !< Restart data unit
  INTEGER, PARAMETER :: dakota_unit = 11  !< Dakota data unit

  !> Counter for the output files
  INTEGER :: output_idx 

  !> Flag to start a run from a previous output:\n
  !> - T     => Restart from a previous output
  !> - F     => Restart from initial condition read from two_phases.inp
  !> .
  LOGICAL :: restart

  ! -- Variables for the namelists LEFT_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcL , u_bcL

  ! -- Variables for the namelists RIGHT_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcR , u_bcR


  NAMELIST / run_parameters / run_name , restart , max_dt ,                     & 
       t_start , t_end , dt_output , solver_scheme ,                            &
       cfl , limiter , theta , reconstr_coeff , n_RK , verbose_level

  NAMELIST / restart_parameters / restart_file

  NAMELIST / newrun_parameters / x0 , xN , comp_cells , riemann_interface

  NAMELIST / left_state / hB_L , u_L

  NAMELIST / right_state / hB_R , u_R

  NAMELIST / left_boundary_conditions / hB_bcL , u_bcL

  NAMELIST / right_boundary_conditions / hB_bcR , u_bcR

  NAMELIST / source_parameters / grav , friction_model

CONTAINS

  !******************************************************************************
  !> \brief Initialization of the variables read from the input file
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE init_param

    USE parameters , ONLY : n_vars

    IMPLICIT none

    LOGICAL :: lexist

    INTEGER :: i

    !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
    run_name = 'default'
    restart = .FALSE.
    max_dt = 1.d-3
    t_start = 0.0
    t_end = 5.0d-4
    dt_output = 1.d-4
    solver_scheme = 'LxF'
    n_RK = 2
    verbose_level = 0


    cfl = 0.5
    limiter(1:n_vars) = 0
    theta=1.0
    reconstr_coeff = 1.0

    !-- Inizialization of the Variables for the namelist restart parameters
    restart_file = ''

    !-- Inizialization of the Variables for the namelist newrun_parameters
    x0 = 0.D0
    xN = 1.D0
    comp_cells = 500
    riemann_interface = 0.5D0

    !-- Inizialization of the Variables for the namelist left_state
    hB_L = 1.D0
    u_L = 0.D0

    !-- Inizialization of the Variables for the namelist right_state
    hB_R = 0.5D0
    u_R = 0.D0

    !-- Inizialization of the Variables for the namelist left boundary conditions

    hB_bcL%flag = 1 
    hB_bcL%value = 0.d0 

    u_bcL%flag = 1 
    u_bcL%value = 0.d0 

    !-- Inizialization of the Variables for the namelist right boundary conditions

    hB_bcR%flag = 1 
    hB_bcR%value = 0.d0 

    u_bcR%flag = 1 
    u_bcR%value = 0.d0 

    !-- Inizialization of the Variables for the namelist source_parameters
    grav = -9.81D0
    friction_model = 'none'

    input_file = 'shallow_water.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       OPEN(input_unit,FILE=input_file,STATUS='NEW')

       WRITE(input_unit, run_parameters )
       WRITE(input_unit, restart_parameters )
       WRITE(input_unit, newrun_parameters )
       WRITE(input_unit, left_state )
       WRITE(input_unit, right_state )
       WRITE(input_unit, left_boundary_conditions )
       WRITE(input_unit, right_boundary_conditions )
       WRITE(input_unit, source_parameters )

       n_batimetry_profile = 2

       ALLOCATE( batimetry_profile(2,n_batimetry_profile) )

       batimetry_profile(1,1) = 0.d0
       batimetry_profile(2,1) = 1.d0

       batimetry_profile(2,2) = 0.d0
       batimetry_profile(2,2) = 10.d0


       WRITE(input_unit,*) '''BATIMETRY_PROFILE'''
       WRITE(input_unit,*) n_batimetry_profile

       DO i = 1, n_batimetry_profile

          WRITE(input_unit,108) batimetry_profile(1:2,i)

108       FORMAT(2(1x,e14.7))


       END DO

       CLOSE(input_unit)

       WRITE(*,*) 'Input file shallow_water.inp not found'
       WRITE(*,*) 'A new one with default values has been created'
       STOP

    ELSE

    END IF

    ! output file index
    output_idx = 0

  END SUBROUTINE init_param

  !******************************************************************************
  !> \brief Read the input file
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_param

    USE parameters, ONLY : bcL , bcR

    IMPLICIT none

    REAL*8 :: max_cfl

    LOGICAL :: tend1 
    CHARACTER(LEN=80) :: card

    INTEGER :: i

    CHARACTER(LEN=4) :: idx_string

    OPEN(input_unit,FILE=input_file,STATUS='old')


    ! ------- READ run_parameters NAMELIST -----------------------------------
    READ(input_unit, run_parameters )

    IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' ) .AND. &
         ( solver_scheme .NE. 'GFORCE' ) ) THEN

       WRITE(*,*) 'WARNING: no correct solver scheme selected',solver_scheme
       WRITE(*,*) 'Chose between: LxF, GFORCE or KT'
       STOP

    END IF

    IF  ( ( solver_scheme.EQ.'LxF' ) .OR. ( solver_scheme.EQ.'GFORCE' ) ) THEN 

       max_cfl = 1.0

    ELSEIF ( solver_scheme .EQ. 'KT' ) THEN

       max_cfl = 0.5

    END IF


    IF ( ( cfl .GT. max_cfl ) .OR. ( cfl .LT. 0.D0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of cfl ',cfl
       WRITE(*,*) 'Choose a value between 0.0 and ',max_cfl
       READ(*,*)

    END IF


    WRITE(*,*) 'Limiters',limiter

    IF ( ( MAXVAL(limiter) .GT. 3 ) .OR. ( MINVAL(limiter) .LT. 0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong limiter ',limiter
       WRITE(*,*) 'Choose among: none, minmod,superbee,van_leer'
       STOP         

    END IF

    IF ( ( reconstr_coeff .GT. 1.0D0 ) .OR. ( reconstr_coeff .LT. 0.D0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of reconstr_coeff ',reconstr_coeff
       WRITE(*,*) 'Change the value between 0.0 and 1.0 in the input file'
       READ(*,*)

    END IF

    IF ( restart ) THEN

       READ(input_unit,restart_parameters)

    ELSE

       READ(input_unit,newrun_parameters)
       READ(input_unit,left_state)
       READ(input_unit,right_state)

    END IF

    READ(input_unit,left_boundary_conditions)

    bcL(1) = hB_bcL 
    bcL(2) = u_bcL 

    READ(input_unit,right_boundary_conditions)

    bcR(1) = hB_bcR 
    bcR(2) = u_bcR 

    READ(input_unit, source_parameters )

    tend1 = .FALSE.

    WRITE(*,*) 'search batimetry_profile'

    batimetry_profile_search: DO

       READ(input_unit,*, END = 200 ) card

       IF( TRIM(card) == 'BATIMETRY_PROFILE' ) THEN

          EXIT batimetry_profile_search

       END IF

    END DO batimetry_profile_search

    READ(input_unit,*) n_batimetry_profile

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'n_batimetry_profile' ,            &
         n_batimetry_profile

    ALLOCATE( batimetry_profile(2,n_batimetry_profile) )

    DO i = 1, n_batimetry_profile

       READ(input_unit,*) batimetry_profile(1:2,i)

       IF ( verbose_level .GE. 1 ) WRITE(*,*) i,batimetry_profile(1:2,i)

    END DO

    GOTO 210
200 tend1 = .TRUE.
210 CONTINUE


    CLOSE( input_unit )

    bak_name = TRIM(run_name)//'.bak'

    OPEN(backup_unit,file=bak_name,status='unknown')

    WRITE(backup_unit, run_parameters )

    IF ( restart ) THEN

       WRITE(backup_unit,restart_parameters)

    ELSE

       WRITE(backup_unit,newrun_parameters)
       WRITE(backup_unit,left_state)
       WRITE(backup_unit,right_state)

    END IF

    WRITE(backup_unit,left_boundary_conditions)
    WRITE(backup_unit,right_boundary_conditions)

    WRITE(backup_unit, source_parameters )

    WRITE(backup_unit,*) '''BATIMETRY_PROFILE'''
    WRITE(backup_unit,*) n_batimetry_profile

    DO i = 1, n_batimetry_profile

       WRITE(backup_unit,107) batimetry_profile(1:2,i)

107    FORMAT(2(1x,e14.7))


    END DO


    CLOSE(backup_unit)


  END SUBROUTINE read_param


  !******************************************************************************
  !> \brief Read the solution from the restart unit
  !
  !> This subroutine is called when the parameter "restart" in the input 
  !> file is TRUE. Then the initial solution is read from a file. 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_solution

    USE geometry, ONLY : init_grid
    USE solver, ONLY : allocate_solver_variables

    USE geometry , ONLY : comp_cells , x0 , dx
    USE parameters, ONLY : n_vars
    USE parameters, ONLY : t_start
    USE solver, ONLY : q

    IMPLICIT none

    INTEGER :: j
    INTEGER :: i

    CHARACTER(LEN=30) :: string

    OPEN(restart_unit,FILE=restart_file,STATUS='old')

    READ(restart_unit,1001) x0,dx,comp_cells,t_start
1001 FORMAT(e18.8,'    x0', /, e18.8,'    dx', /, i18,'    cells', /, e18.8,'    t', /)


    CALL init_grid

    CALL allocate_solver_variables


    DO i = 1,n_vars
       DO j = 1,comp_cells

          ! Exponents with more than 2 digits cause problems reading
          ! into matlab... reset tiny values to zero:
          IF ( dabs(q(i,j)) .LT. 1d-99) q(i,j) = 0.d0

       END DO

       READ(restart_unit,1003) (q(i,j), j=1,n_vars)
1003   FORMAT(4e20.12)

    END DO

    j = SCAN(restart_file, '.' , .TRUE. )

    string = TRIM(restart_file(j+2:j+5))

    READ( string,* ) output_idx

  END SUBROUTINE read_solution

  !******************************************************************************
  !> \brief Write the solution on the output unit
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !
  !> \param[in]   t      output time
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE output_solution(t)

    ! external procedures
    USE constitutive, ONLY : qc_to_qp

    ! external variables
    USE constitutive, ONLY : h , u
    USE geometry , ONLY : comp_cells , x0 , dx , B_cent
    USE parameters, ONLY : n_vars
    USE parameters, ONLY : t_output , dt_output
    USE solver, ONLY : q

    IMPLICIT none

    REAL*8, INTENT(IN) :: t

    CHARACTER(LEN=4) :: idx_string

    REAL*8 :: qp(n_vars)

    INTEGER :: j
    INTEGER :: i

    output_idx = output_idx + 1

    idx_string = lettera(output_idx-1)

    output_file = TRIM(run_name)//'.q'//idx_string

    WRITE(*,*) 'WRITING ',output_file

    OPEN(output_unit,FILE=output_file,status='unknown',form='formatted')

    WRITE(output_unit,1002) x0,dx,comp_cells,t

    DO j = 1,comp_cells
       DO i = 1,n_vars

          ! Exponents with more than 2 digits cause problems reading
          ! into matlab... reset tiny values to zero:
          IF ( dabs(q(i,j)) .LT. 1d-99) q(i,j) = 0.d0

       END DO

       WRITE(output_unit,1005) (q(i,j), i=1,n_vars)

    END DO

    WRITE(output_unit,*) ' '
    WRITE(output_unit,*) ' '

    CLOSE(output_unit)

    output_file = TRIM(run_name)//'.p'//idx_string

    WRITE(*,*) 'WRITING ',output_file

    OPEN(output_unit,FILE=output_file,status='unknown',form='formatted')

    WRITE(output_unit,1002) x0,dx,comp_cells,t

    DO j = 1,comp_cells

       DO i = 1,n_vars

          ! Exponents with more than 2 digits cause problems reading
          ! into matlab... reset tiny values to zero:
          IF ( dabs(q(i,j)) .LT. 1d-99) q(i,j) = 0.d0

       END DO

       CALL qc_to_qp(q(:,j),B_cent(j),qp(:))

       WRITE(output_unit,1005) REAL(h) , REAL(u) , B_cent(j) ,                &
            REAL(h) + B_cent(j)

    END DO

    WRITE(output_unit,*) ' '
    WRITE(output_unit,*) ' '

    CLOSE(output_unit)

1002 FORMAT(e18.8,'    x0', /, e18.8,'    dx', /, i18,'    cells', /, e18.8,'    t', /)


1005 FORMAT(4e20.12)


    t_output = t + dt_output

  END SUBROUTINE output_solution


  SUBROUTINE close_units

    IMPLICIT NONE

  END SUBROUTINE close_units

  !******************************************************************************
  !> \brief Numeric to String conversion
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !
  !> \param[in]   k      integer to convert         
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  CHARACTER*4 FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera

  SUBROUTINE output_dakota

    IMPLICIT NONE

    dakota_file = 'shallow_water.out'

    OPEN(dakota_unit,FILE=dakota_file,status='unknown',form='formatted')

    CLOSE(dakota_unit)

  END SUBROUTINE output_dakota


END MODULE inpout

