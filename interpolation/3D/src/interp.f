      module interp_mod 
      implicit none
      private

      public :: interpd
      public :: neville, bounds, boundsN, xyinterp

      contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       interpd -- general driver for interpolation in 3D
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        real*8 function interpd(Nx, Ny, Nz, dx, dy, dz, xmin, xmax,
     &                          ymin, ymax, zmin, zmax, x0, y0, z0, f)

        integer Nx, Ny, Nz
        real*8 dx, dy, dz
        real*8 xmin, xmax, ymin, ymax, zmin, zmax
        real*8 x0, y0, z0
        real*8, dimension(Nx,Ny,Nz) :: f

        integer i
        integer zlowerN, zupperN
        real*8, dimension(3) :: uppers
        integer, dimension(3) :: uppersN
        real*8, dimension(:), allocatable :: zvals, fvals

        uppers = bounds(dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax,
     &                  x0, y0, z0)
        uppersN = boundsN(dx, dy, dz, xmin, ymin, zmin, uppers)

        zlowerN = uppersN(3)-4
        zupperN = uppersN(3)

        allocate(zvals(zlowerN:zupperN))
        allocate(fvals(zlowerN:zupperN))

        do i = zupperN, zlowerN, -1
          zvals(i) = zmin + (i-1)*dz
          fvals(i) = xyinterp(dx, dy, xmin, ymin, Nx, Ny, Nz, uppersN,
     &                        f, i, x0, y0)
        end do

        interpd = neville(zvals,fvals,z0)

        return
        end function interpd

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       neville -- interpolates a value at a point via the Neville
c                  algorithm. Note that the interpolation is hard-coded
c                  to be fourth-order. Input arguments are a list of 5
c                  axis values, a list of 5 known function values at
c                  those points, and a point at which the interpolation
c                  is desired
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        real*8 function neville(x, f, p)

        implicit none

        real*8, intent(in), dimension(0:4) :: x
        real*8, intent(in), dimension(0:4) :: f
        real*8, intent(in) :: p

        real*8, dimension(0:4,0:4) :: Q = 0d0
        integer i, j

        do i=0,4
          Q(i,0) = f(i)
        end do

        ! Neville's algorithm to fourth-order
        do i=1,4
        do j=1,i
          Q(i,j)=(p-x(i-j))*Q(i,j-1)-(p-x(i))*Q(i-1,j-1)
          Q(i,j)=Q(i,j)/(x(i)-x(i-j))

          ! prevent underflow
          if ( Q(i,j) .lt. 1d-50 ) then
            neville = 0d0
          end if
        end do
        end do

c        write(*,*) x(0), Q(0,0)
c        write(*,*) x(1), Q(1,0), Q(1,1)
c        write(*,*) x(2), Q(2,0), Q(2,1), Q(2,2)
c        write(*,*) x(3), Q(3,0), Q(3,1), Q(3,2), Q(3,3)
c        write(*,*) x(4), Q(4,0), Q(4,1), Q(4,2), Q(4,3), Q(4,4)

        neville = Q(4,4)

        return
        end function neville

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       bounds -- returns the upper (x,y,z) bounds of the relevant
c                 interpolating cube
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function bounds(dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax,
     &                  x0, y0, z0)
     &  result(uppers) 
        implicit none
  
        real*8, intent(in) :: dx, dy, dz
        real*8, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
        real*8, intent(in) :: x0, y0, z0

        real*8, dimension(3) :: uppers

        real*8 xupper, yupper, zupper
  
        if ( x0+4d0*dx .GE. xmax ) then
          xupper = xmax
        else
          xupper = xmin + real(int((x0-xmin)/dx))*dx + 4d0*dx
        end if
  
        if ( y0+4d0*dy .GE. ymax ) then
          yupper = ymax
        else
          yupper = ymin + real(int((y0-ymin)/dy))*dy + 4d0*dy
        end if
  
        if ( z0+4d0*dz .GE. zmax ) then
          zupper = zmax
        else
          zupper = zmin + real(int((z0-zmin)/dz))*dz + 4d0*dz
        end if

        uppers = [xupper, yupper, zupper]

        return
        end function bounds

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       boundsN -- returns the grid point number for an (x,y,z) tuple
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function boundsN(dx, dy, dz, xmin, ymin, zmin, uppers)
     &  result(uppersN) 
        implicit none
  
        real*8, intent(in) :: dx, dy, dz
        real*8, intent(in) :: xmin, ymin, zmin
        real*8, intent(in), dimension(3) :: uppers

        integer, dimension(3) :: uppersN

        integer xupperN, yupperN, zupperN
  
        xupperN = idnint((uppers(1)-xmin)/dx) + 1
        yupperN = idnint((uppers(2)-ymin)/dy) + 1
        zupperN = idnint((uppers(3)-zmin)/dz) + 1
  
        uppersN = [xupperN, yupperN, zupperN]

        return
        end function boundsN

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       xyinterp -- interpolates along lines of known x then uses those
c                   values to interpolate along y and return the
c                   interpolation estimate of f(x0,y0) in the 2D plane
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function xyinterp(dx, dy, xmin, ymin, Nx, Ny, Nz, uppersN, f,
     &                    zindex, x0, y0)
     &  result(xyresult) 
        implicit none
  
        real*8, intent(in) :: dx, dy
        real*8, intent(in) :: xmin, ymin
        integer, intent(in) :: Nx, Ny, Nz
        integer, intent(in), dimension(3) :: uppersN 
        real*8, intent(in), dimension(Nx,Ny,Nz) :: f
        integer, intent(in) :: zindex
        real*8, intent(in) :: x0, y0

        real*8 xyresult

        integer i, j, k, m
        real*8, dimension(5) :: xresult
        real*8, dimension(5) :: x, y, g

        m = 1
        do j = uppersN(2)-4, uppersN(2)
        y(m) = ymin + (j-1)*dy
        k = 1
        do i = uppersN(1)-4, uppersN(1)
          x(k) = xmin + (i-1)*dx
          g(k) = f(i,j,zindex)
          k = k + 1
        end do
        xresult(m) = neville(x,g,x0)
        m = m + 1
        end do

        xyresult = neville(y,xresult,y0)

        return
        end function xyinterp

      end module interp_mod
