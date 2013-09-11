!=============================================================================!
! This subroutine replaces the DUMSES init subroutine.                        !
!                                                                             !
! The DUMSES init subroutine takes 5 arguments.  The first three are          !
! parameters for DUMSES output (e.g., the number for the current output       !
! file).  These are unnecessary because we are not using the DUMSES output.   !
! The next two are the current processor ID and the number of processors.     !
! Since we are writing this program to be serial, there is only one processor !
! and it needs no ID.                                                         !
!=============================================================================!

subroutine init()

   use hydro_parameters
   use variables

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! ==========================================================================
   ! Initialize

   call read_input()

   ! Calculate state vector sizes
   iu1 = 1 - nghost
   iu2 = nx + nghost

   ! Calculate flux vector sizes
   if1 = 1
   if2 = nx + 1

   ! Allocate arrays
   allocate(uin         (iu1:iu2,nvar+3)); uin          = 0.0d0
   allocate(old_uin     (iu1:iu2,nvar+3)); old_uin      = 0.0d0
   allocate(gravin      (iu1:iu2       )); gravin       = 0.0d0
   allocate(flux        (iu1:iu2,nvar  )); flux         = 0.0d0
   allocate(flux_pre_tot(iu1:iu2,nvar  )); flux_pre_tot = 0.0d0
   !allocate(emfx        (iu1:iu2,nvar  )); emfx         = 0.0d0
   !allocate(emfy        (iu1:iu2,nvar  )); emfy         = 0.0d0
   !allocate(emfz        (iu1:iu2,nvar  )); emfz         = 0.0d0

   ! Times
   dt   = 0.0d0

   ! Create grid
   allocate(x(iu1:iu2))
   dx = (xmax - xmin) / nx
   x = xmin + 0.5*dx + (/ (i-1,i=iu1,iu2) /)*dx

   ! Initialize the problem
   call condinit()
   call bval(uin,dt,dx,bval_type)

   ! Initialize geometry (cell volumes and face areas)
   call init_geom()

end subroutine init

!=============================================================================!
! This subroutine replaces the DUMSES read_input subroutine.                  !
!                                                                             !
!=============================================================================!

subroutine read_input()

   use hydro_parameters
   use variables

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer, parameter :: LUN = 17

   ! Namelists
   namelist /scheme_params/nx_glob,riemann,slope_type,bval_in,bval_out,&
                           &courant,fargo
   namelist /model_params/geom,xmin,xmax,gamma,nu,eta

   ! ==========================================================================
   ! Set the default parameters

   ! scheme_params
   nx_glob    = 1
   riemann    = "hlld"
   slope_type = 1
   bval_in    = 1
   bval_out   = 1
   courant    = 0.8
   fargo      = .false.

   ! model_params
   geom   = 1
   xmin   = 0.0d0
   xmax   = 1.0d0
   gamma  = 1.001d0
   nu     = 0.0d0
   eta    = 0.0d0

   ! ==========================================================================
   ! Read the input file

   open(unit=LUN, file="input", status="old")
   ! Don't need start params
   open(LUN, scheme_params)
   open(LUN, model_params)
   ! Don't need output params
   close(LUN)

   ! ==========================================================================
   ! Some repackaging and processing

   nx = nx_glob

   bval_type(1) = bval_in
   bval_type(2) = bval_out

   cartesian   = .false.
   cylindrical = .false.
   spherical   = .false.
   select case (geom)
   case (1)
      cartesian = .true.
   case (2)
      cylindrical = .true.
   case (3)
      spherical = .true.
   end select

   ! ==========================================================================
   ! Solver details

   select case (riemann)
   case ('roe')
      iriemann = iroe
   case ('llf')
      iriemann = illf
   case ('hll')
      iriemann = ihll
   case ('hlld')
      iriemann = ihlld
   case ('upwind')
      iriemann = iupwind
   case ('hydro')
      iriemann = iacoustic
   end select

end subroutine read_input

!=============================================================================!
! This subroutine replaces the DUMSES init_geom subroutine.                   !
!                                                                             !
!=============================================================================!

subroutine init_geom()

   use hydro_parameters
   use variables

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer :: ilo, ihi

   ! ==========================================================================
   ! Set up geometrical factors

   allocate(dv(iu1:iu2)); dv = 0.0d0
   allocate(ds(iu1:iu2)); ds = 0.0d0

   ilo = min(1, iu1+1)
   ihi = max(1, iu2-1)

   if cartesian then
      do i = ilo, ihi
         dv(i) = dx
         ds(i) = 1
      end do
   else if cylindrical then
      do i = ilo, ihi
         dv(i) = dx * x(i)
         ds(i) = 0.5d0*(x(i)+x(i-1))
      end do
   else if spherical then
      do i = ilo, ihi
         dv(i) = dx * x(i) * x(i)
         ds(i) = 0.5d0*(x(i)+x(i-1))
      end do
   end if

end subroutine init_geom

!=============================================================================!
! This subroutine allows the user to do some initialization.                  !
!                                                                             !
!=============================================================================!

subroutine user_init()

   implicit none

   ! ==========================================================================
   ! Declare variables

end subroutine user_init
