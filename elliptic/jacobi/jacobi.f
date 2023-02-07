c=======================================================================
c
c     jacobi: solves a linear system using Jacobi iteration
c
c     The problem solved is the Laplace equation on [0,1]x[0,1]:
c
c           u(x,0)=0, u(x,1)=f(x),    0 < x < 1
c           u(0,y)=0, u(1,y)=0,       0 < y < 1
c
c     where
c                  { 75x,       if x<= 2/3
c           f(x) = {
c                  { 150(1-x),  if x > 2/3
c
c     and the exact solution is 
c
c                   450  infty / sin(2pi*n/3)sin(n*pi*x)*sinh(n*pi*y) \
c           u(x,y)= ---   SUM |  -----------------------------------   |
c                   pi^2  n=1  \           n^2 * sinh(n*pi)           /
c
c=======================================================================

      program jacobi

      implicit none

      real*8        pi
      parameter    (pi=3.14159265358979)

c     2D allocatable arrays to store the grid values and source
      real*8, allocatable, dimension(:,:) ::    u,  unew, f

      integer       i,  j,  k     ! loop variables
      real*8        x,  y,  m     ! for converting loop vars to reals

      integer       n             ! grid size (n x n)
      real*8        h             ! discretization scale
      integer       iter          ! current Jacobi iteration
      integer       maxiter       ! max # of Jacobi iterations
      real*8        err           ! l-infinity norm
      real*8        tol           ! error tolerance (stopping criteria)
      real*8        abserr, exact ! for absolute error of numerical sol

      integer       unum, stat    ! used to read/open/close a file
      integer       unum1, unum2  ! used to read/open/close a file
      integer       getu          ! defined in p410f library

c      write(*,*) "Enter the desired size of the grid: "
c      read(*,*) n
      n=100  ! grid will be (n+1) by (n+1)

c      write(*,*) "Enter the desired grid spacing: "
c      read(*,*) h
      h=0.01

c     Allocate the arrays (start indexing at 0 for simplicity)
      allocate( u(0:n,0:n) )
      allocate( unew(0:n,0:n) )
      allocate( f(0:n,0:n) )

c     Parameters for the algorithm
      maxiter=50000
      tol=0.00000001    
     
c     Set the initial guess and boundary conditions
      do i=0,n
      do j=0,n
        u(i,j)=0.0     ! initial guess is all zeros
        f(i,j)=0.0     ! no source function

        x=real(i)*h    ! y=j*h, x=i*h
        if (j .EQ. n ) then 
          if ( x .LE. 0.666) then
            u(i,j)=75.0*x
          else
            u(i,j)=150.0*(1.0-x)
          endif
        endif

      end do
      end do
      unew=u

      write(*,*)
      write(*,*) 'Maximum iterations: ', maxiter
      write(*,*) 'Specified tolerance:', tol 
      write(*,*)

      do iter=1,maxiter
         err = 0d0   
         do i=1,n-1
           do j=1,n-1   ! using 5-point stencil:
              unew(i,j)=0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)
     &                  - h**2*f(i,j));
           end do
         end do
         do i=0,n
         do j=0,n
           if ( abs(unew(i,j)-u(i,j)) > err ) err=abs(unew(i,j)-u(i,j))
         end do
         end do
         if ( err < tol ) exit
         u=unew   
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
      unum1 = getu()  ! for writing exact solution
      unum2 = unum1+1 ! for writing absolute error
      open(unit=unum1, file='solution_exact', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum2, file='solution_abserr', status='replace',
     &     action='write',iostat=stat)
      do i=0,n
      do j=0,n
        x=real(i)*h  ! x=i*h
        y=real(j)*h  ! y=j*h
        exact=0.0d0
        do k=1,200  ! summation (n=1...infty) in exact solution 
          m=real(k)
          exact=exact+(sin(2.0*m*pi/3.0)*sin(m*pi*x)
     &               *sinh(m*pi*y))/(m**2*sinh(m*pi))
        end do
        exact=(450.0/pi**2.0)*exact
        abserr=exact
        abserr=abs(u(i,j)-exact)
        write(unum1,400,iostat=stat) x, y, exact
        write(unum2,400,iostat=stat) x, y, abserr
      end do
      end do
      close(unit=unum1, status='keep', iostat=stat)
      close(unit=unum2, status='keep', iostat=stat)
      
      stop
      end
