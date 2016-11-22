!********************************************************************************
!> \mainpage   IMEX_SfloW - Shallow Water Finite volume solver
!> IMEX_SfloW is a FORTRAN90 code designed to solve an hyperbolic 
!> system of partial differential equations with relaxation and source
!> terms. 
!> The model is discretized in time with an explicit-implicit Runge-Kutta
!> method where the hyperbolic part is solved explicetely and the other
!> terms (relaxation and surce) are treated implicitely.\n
!> The finite volume solver for the hyperbolic part of the system is based
!> on a semidiscrete central scheme and it is not tied on the specific 
!> eigenstructure of the model.\n
!> The implicit part is solved with a Newton-Raphson method where the 
!> elements of the Jacobian of the nonlinear system are evaluated 
!> numerically with a complex step derivative technique.
!> Version 1.0:\n
! 
!> Github project page: http://demichie.github.io/IMEX_SfloW/
!> \n
!> \authors Mattia de' Michieli Vitturi (*,**)
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!> (**) School of Earth and Space Exploration, Arizona State University \n
!>      Tempe, AZ, USA \n
!>      E-mail: mdemichi@asu.edu \n
!********************************************************************************

!> \brief Main Program 
PROGRAM IMEX_SfloW

  USE constitutive, ONLY : init_problem_param

  USE geometry, ONLY : init_grid
  USE geometry, ONLY : B_cent , dx

  USE init, ONLY : riemann_problem, initial_conditions

  USE inpout, ONLY : init_param
  USE inpout, ONLY : read_param
  USE inpout, ONLY : output_solution
  USE inpout, ONLY : read_solution
  USE inpout, ONLY : close_units
  USE inpout, ONLY : output_dakota

  USE solver, ONLY : allocate_solver_variables
  USE solver, ONLY : deallocate_solver_variables
  USE solver, ONLY : imex_RK_solver
  USE solver, ONLY : timestep

  USE inpout, ONLY : restart

  USE parameters, ONLY : t_start
  USE parameters, ONLY : t_end
  USE parameters, ONLY : t_output
  USE parameters, ONLY : riemann_flag

  USE solver, ONLY : q , dt

  IMPLICIT NONE

  REAL*8 :: t
  REAL*8 :: t1 , t2

  CALL cpu_time(t1)

  CALL init_param

  CALL read_param

  IF ( restart ) THEN

     CALL read_solution

  ELSE

     CALL init_grid

     CALL init_problem_param

     CALL allocate_solver_variables

     ! riemann problem defined in file.inp
     IF(riemann_flag.EQV..TRUE.)THEN

        CALL riemann_problem

     ! generic problem defined by initial conditions function (in init.f90)
     ELSE

        CALL initial_conditions

     ENDIF

  END IF

  t = t_start

  WRITE(*,*) 't_start =',t

  CALL output_solution(t)

  DO WHILE ( t .LT. t_end )

     CALL timestep

     IF ( t+dt .GT. t_end ) dt = t_end - t
     IF ( t+dt .GT. t_output ) dt = t_output - t

     CALL imex_RK_solver

     t = t+dt

     WRITE(*,*) 't =',t,' dt =',dt,' h_tot =',dx*(SUM(q(1,:)-B_cent(:)))

     IF ( ( t .GE. t_output ) .OR. ( t .GE. t_end ) ) CALL output_solution(t)

  END DO

  CALL output_dakota

  CALL deallocate_solver_variables

  CALL close_units

  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'

END PROGRAM IMEX_SFLOW

