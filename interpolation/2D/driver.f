cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MODULES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module globals
        implicit none
  
        real*8 pi
        parameter ( pi = 3.14159265358979d0 )
      end module


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MAIN PROGRAM
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program driver
      implicit none

      integer :: i, j
      real*8  :: ansy, ansy1, ansy2
      real*8  :: x1, x1l, x1u, x2, x2l, x2u
      real*8, dimension(4) :: y, y1, y2, y12

      integer :: Nx1, Nx2
      real*8  :: dx1, dx2
      real*8  :: x1min, x1max, x2min, x2max
      real*8  :: f, f_x, f_y, f_xy

      integer :: iscale

      ! grid properties
      x1min = 0d0
      x1max = 1d0 
      x2min = 0d0
      x2max = 1d0
      Nx1 = 16
      Nx2 = 16
      dx1 = (x1max-x1min)/dble(Nx1)
      dx2 = (x2max-x2min)/dble(Nx2)

      ! how many times to interpolate between grid points
      iscale = 16 

      do i=0,Nx1*iscale
      do j=0,Nx2*iscale

        ! compute desired value
        x1 = i*dx1/iscale
        x2 = j*dx2/iscale

        ! get closest grid points to desired value
        x1l = dx1*floor(x1/dx1)
        x1u = x1l+dx1
        x2l = dx2*floor(x2/dx2)
        x2u = x2l+dx2
        if (x1 .GE. x1max) then 
          x1l = x1l-dx1
          x1u = x1u-dx1
        endif
        if (x2 .GE. x2max) then 
          x2l = x2l-dx2
          x2u = x2u-dx2
        endif

        write(*,*)
        write(*,'(A,ES22.15,ES22.15)') "Desired point:  ", x1,  x2
        write(*,'(A,ES22.15,ES22.15)') "x1l, x1u:       ", x1l, x1u 
        write(*,'(A,ES22.15,ES22.15)') "x2l, x2u:       ", x2l, x2u 
        
        ! set function and derivative values
        y(1)   = f(x1l,x2l)
        y(2)   = f(x1u,x2l)
        y(3)   = f(x1u,x2u)
        y(4)   = f(x1l,x2u)
        y1(1)  = f_x(x1l,x2l)
        y1(2)  = f_x(x1u,x2l)
        y1(3)  = f_x(x1u,x2u)
        y1(4)  = f_x(x1l,x2u)
        y2(1)  = f_y(x1l,x2l)
        y2(2)  = f_y(x1u,x2l)
        y2(3)  = f_y(x1u,x2u)
        y2(4)  = f_y(x1l,x2u)
        y12(1) = f_xy(x1l,x2l)
        y12(2) = f_xy(x1u,x2l)
        y12(3) = f_xy(x1u,x2u)
        y12(4) = f_xy(x1l,x2u)

        write(*,'(A,ES22.15)') "f(x1l, x2l):    ", y(1)
        write(*,'(A,ES22.15)') "f(x1l, x2u):    ", y(2)
        write(*,'(A,ES22.15)') "f(x1u, x2u):    ", y(3)
        write(*,'(A,ES22.15)') "f(x1u, x2l):    ", y(4)

        ! interpolate
        call interp(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
 
        write(*,'(A,ES22.15)') "exact value:        ", f(x1,x2)
        write(*,'(A,ES22.15)') "interpolated value: ", ansy
        write(*,'(A,ES22.15)') "error:              ", f(x1,x2)-ansy
        write(*,'(A,ES22.15)') 

        ! write interpolated solution to file
        open(unit=20, file="interp.dat", access="append",
     &       status="unknown")
          write(20, *) x1, x2, ansy
        close(20)

        ! write exact solution to file
        open(unit=20, file="exact.dat", access="append",
     &       status="unknown")
          write(20, *) x1, x2, f(x1,x2)
        close(20)

        ! write error to file
        open(unit=20, file="error.dat", access="append",
     &       status="unknown")
          write(20, *) x1, x2, f(x1, x2)-ansy
        close(20)

      end do
      end do

      end program driver


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function f(x,y)
        use globals
          real*8 :: f
          real*8 :: x, y
          f = sin(2d0*x*pi)*cos(2d0*y*pi)
          return
      end function

      function f_x(x,y)
        use globals
          real*8 :: f_x
          real*8 :: x, y
          f_x = 2d0*pi*cos(2d0*x*pi)*cos(2d0*y*pi)
          return
      end function

      function f_y(x,y)
        use globals
          real*8 :: f_y
          real*8 :: x, y
          f_y = -2d0*sin(2d0*x*pi)*pi*sin(2d0*y*pi)
          return
      end function

      function f_xy(x,y)
        use globals
          real*8 :: f_xy
          real*8 :: x, y
          f_xy = -4d0*pi**2d0*cos(2d0*x*pi)*sin(2d0*y*pi)
          return
      end function
