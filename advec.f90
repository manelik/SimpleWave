!
!
!
  program advec

! This is a simple code for the advection equation

! define:
!
! pi  = dphi/dt
! psi = dphi/dx


! *************************
! ***   WAVE EQUATION   ***
! *************************

! Declare variables.

  implicit none

  integer i,j,l         ! Counters
  integer Nt            ! Total number of time steps.
  integer Nx            ! Total number of grid points.
  integer Noutput       ! Frequency of output.

  real(8) dx            ! Grid spacing.
  real(8) dt            ! Time step.
  real(8) t             ! Time.

  real(8) rho       ! courant parameter

  real(8) v             ! Wave speed.

  real(8) x0            ! Center of initial gaussian.
  real(8) s0            ! Width of initial gaussian.
  real(8) a0            ! Amplitude of initial gaussian.

  character(20) method  ! Integration method.

  real(8), allocatable, dimension(:) :: x       ! Position.

  real(8), allocatable, dimension(:) :: phi     ! Wave function.
  real(8), allocatable, dimension(:) :: phi_p   ! Old phi.
  real(8), allocatable, dimension(:) :: sphi    ! Source for phi.



! **************************
! ***   GET PARAMETERS   ***
! **************************

  print *
  print *, 'Give grid spacing dx'
  read(*,*) dx

  print *
  print *, 'Give time step dt'
  read(*,*) dt

  print *
  print *, 'Give wave speed'
  read(*,*) v

  print *
  print *, 'Give total number of grid points Nx'
  read(*,*) Nx

  print *
  print *, 'Give total number of time steps'
  read(*,*) Nt

  print *
  print *, 'Give frequency of output'
  read(*,*) Noutput

  print *
  print *, 'Give integration method (only f_euler for now)'
  read(*,*) method


! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

  allocate(x(0:Nx))

  allocate(phi(0:Nx),phi_p(0:Nx),sphi(0:Nx))


! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

  do i=0,Nx
     x(i) = dble(i)*dx
  end do


! ****************************
! ***   OUTPUT TO SCREEN   ***
! ****************************

  print *,'------------------------------'
  print *,'|  Time step  |     Time     |'
  print *,'------------------------------'


! ************************
! ***   INITIAL DATA   ***
! ************************

! Initialize time.

  t = 0.0D0

! Initialize courant factor

  rho = dt/dx

! Parameters for initial data.

  a0 = 1.0D0
  x0 = dble(Nx/2)*dx
  s0 = 1.0D0

! Initial data (gaussian).

  do i=0,Nx
     phi(i) = a0*exp(-(x(i)-x0)**2/s0**2)
  end do

! Output to screen.

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *****************************
! ***   OPEN OUTPUT FILES   ***
! *****************************

   open(1,file='phi.xl',form='formatted',status='replace')


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************

! Wave function.

  write(1,"(A8,ES14.6)") '"Time = ',t

  do i=0,Nx
     if (dabs(phi(i)) > 1.0D-50) then
        write(1,"(2ES16.8)") x(i),phi(i)
     else
        write(1,"(2ES16.8)") x(i),0.0D0
     end if
  end do


! Leave blank spaces before next time level.

  write(1,*)
  write(1,*)



! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

     phi_p = phi

!    Forward Euler

     if (method=='f_euler') then 

        do j = 1, (Nx-1)
           phi(j) = phi_p(j) - 0.5D0*v*rho*( phi_p(j+1)-phi_p(j-1) )
        end do

!    We leave boundaries alone for the moment

!    Unknown.

     else

        print *, 'Unknown integartion method.'
        print *
        stop

     end if


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

     if (mod(l,Noutput).eq.0) then

!       Wave function.

        write(1,"(A8,ES14.6)") '"Time = ',t

        do i=0,Nx
           if (dabs(phi(i)) > 1.0D-50) then
              write(1,"(2ES16.8)") x(i),phi(i)
           else
              write(1,"(2ES16.8)") x(i),0.0D0
           end if 
        end do


!       Leave blank spaces before next time level.

        write(1,*)
        write(1,*)

     end if


!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

     if (mod(l,Noutput).eq.0) then
        write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',l,'   | ',t,'  | '
     end if

  end do

  print *,'------------------------------'


! ******************************
! ***   CLOSE OUTPUT FILES   ***
! ******************************

  close(1)


! ***************
! ***   END   ***
! ***************

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *
  print *

  end program advec

