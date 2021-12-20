c=======================================================================
c
c     sor: solves the Poisson equation using successive overrelaxation
c          in Cartesian coordinates
c
c     The problem solved is the Poisson equation on [0,1]x[0,1]:
c
c           del u(x,y) = f(x,y)
c           u(x,0) = 0, u(x,1) = 0,      0 < x < 1
c           u(0,y) = 0, u(1,y) = y(1-y), 0 < y < 1
c
c       where
c
c           f(x,y) = 6xy(1-y)-2x^3
c
c       and the exact solution is 
c
c           u(x,y) = y(1-y)x^3
c
c=======================================================================

      program sor

      implicit none

      real*8        pi
      parameter    (pi=3.14159265358979)

c     2D allocatable arrays to store the grid values and source
      real*8, allocatable, dimension(:,:) ::    u,  uold, f

      integer       i,  j,  k     ! loop variables
      real*8        x,  y,  m     ! for converting loop vars to reals

      integer       n             ! grid size (n x n)
      real*8        h             ! discretization scale
      real*8        w             ! overrelaxation parameter
      integer       iter          ! current iteration
      integer       maxiter       ! max # of iterations
      real*8        err           ! l-infinity norm
      real*8        tol           ! error tolerance (stopping criteria)
      real*8        abserr, exact ! for absolute error of numerical sol

      integer       unum, stat    ! used to read/open/close a file
      integer       getu          ! defined in p410f library
      integer       unum1, unum2, unum3 ! for storing unit numbers

c      write(*,*) "Enter the desired size of the grid: "
c      read(*,*) n
      n=100  ! grid will be (n+1) by (n+1) to cover the entire domain

c      write(*,*) "Enter the desired grid spacing: "
c      read(*,*) h
      h=0.01

c     Allocate the arrays (start indexing at 0 for simplicity)
      allocate( u(0:n,0:n) )
      allocate( uold(0:n,0:n) )
      allocate( f(0:n,0:n) )

c     Parameters for the algorithm
      maxiter=1000
      tol=1e-5  
      w=0.5*(1.0-pi/(n+1))  ! Haberman Exercise 6.6.2
      
      write(*,*)
      write(*,*) 'Maximum iterations:  ', maxiter
      write(*,*) 'Specified tolerance:      ', tol 
      write(*,*) 'Overrelaxation parameter: ', w 
      write(*,*)

c     Set the initial guess
      do i=0,n
      do j=0,n
        x=real(i)*h    ! x=i*h
        y=real(j)*h    ! y=j*h
        u(i,j)=0.0     ! initial guess is all zeros
        f(i,j)=6.0*x*y*(1.0-y)-2.0*x**3.0
      end do
      end do

      do iter=1,maxiter
         uold=u
         err = 0d0   

c        boundary condition at x=1
         do j=0,n
           y=real(j)*h  ! y=j*h
           u(n,j)=y*(1.0-y)
         end do

c        interior points
         do i=1,n-1
           do j=1,n-1   ! using 5-point stencil:
             u(i,j)=u(i,j)+w*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)
     &              -4.0*u(i,j)-h**2*f(i,j))
           end do
         end do
         do i=0,n
         do j=0,n
           if ( abs(uold(i,j)-u(i,j)) > err ) err=abs(uold(i,j)-u(i,j))
         end do
         end do
         if ( err < tol ) exit
      end do

      if ( iter < maxiter ) then
        write(*,200) iter
 200    format( 'Tolerance satisfied in ', I6, ' iterations.' )
      else
        write(*,300) err
 300    format( 'Maxiter reached. Program exited with final error',
     &           F10.5 )
      endif

      write(*,*)

c     Open, write and close the file for the solution
      unum = getu()
      open(unit=unum, file='solution_numerical', status='replace',
     &     action='write',iostat=stat)
      do i=0,n
      do j=0,n
        write(unum,400,iostat=stat) i*h, j*h, u(i,j)
      end do
      end do
 400  format( F5.3, "  ", F5.3, "  ", F12.7 )
      close(unit=unum, status='keep', iostat=stat)

c     Calculate the exact solution and absolute error and write to file
c     Also write the source function to file
      unum1 = getu()  ! for writing exact solution
      unum2 = unum1+1 ! for writing absolute error
      unum3 = unum1+2 ! for writing the source
      open(unit=unum1, file='solution_exact', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum2, file='solution_abserr', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum3, file='solution_source', status='replace',
     &     action='write',iostat=stat)
      do i=0,n
      do j=0,n
        x=real(i)*h  ! x=i*h
        y=real(j)*h  ! y=j*h
        exact=y*(1.0-y)*x**3.0
        abserr=abs(u(i,j)-exact)
        write(unum1,400,iostat=stat) x, y, exact
        write(unum2,400,iostat=stat) x, y, abserr
        write(unum3,400,iostat=stat) x, y, f(i,j)
      end do
      end do
      close(unit=unum1, status='keep', iostat=stat)
      close(unit=unum2, status='keep', iostat=stat)
      close(unit=unum3, status='keep', iostat=stat)

c     Now in 1D, output the exact solution, numerical, and error
      unum1 = getu()  ! for writing exact solution
      unum2 = unum1+1 ! for writing numerical solution
      unum3 = unum1+2 ! for writing absolute error
      open(unit=unum1, file='solution_exact_1D', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum2, file='solution_numerical_1D', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum3, file='solution_abserr_1D', status='replace',
     &     action='write',iostat=stat)
      j=n/2 
      y=real(j)*h  ! y=j*h
      do i=0,n
        x=real(i)*h  ! x=i*h
        exact=y*(1.0-y)*x**3.0
        abserr=abs(u(i,j)-exact)
        write(unum1,400,iostat=stat) x, exact
        write(unum2,400,iostat=stat) x, u(i,j)
        write(unum3,400,iostat=stat) x, abserr
      end do
      close(unit=unum1, status='keep', iostat=stat)
      close(unit=unum2, status='keep', iostat=stat)
      close(unit=unum3, status='keep', iostat=stat)

      stop
      end
