      program driver
      use interp_mod

      implicit none

      integer Nx, Ny, Nz
      real*8 dx, dy, dz
      real*8 xmin, xmax, ymin, ymax, zmin, zmax
      real*8 x0, y0, z0
      real*8 x, y, z
      integer i, j, k
      real*8, dimension(:,:,:), allocatable :: f
      real*8 interpval

      Nx = 1025
      Ny = 1025
      Nz = 1025
      xmin = -4d0
      xmax = 8d0
      ymin = 2d0
      ymax = 8d0
      zmin = 0d0
      zmax = 8d0
      dx = (xmax-xmin)/real(Nx-1)
      dy = (ymax-ymin)/real(Ny-1)
      dz = (zmax-zmin)/real(Nz-1)

      ! construct a 3D surface to interpolate
      allocate(f(Nx, Ny, Nz))
      do i = 1, Nx
        do j = 1, Ny
          do k = 1, Nz
            x = xmin + real(i - 1) * dx
            y = ymin + real(j - 1) * dy
            z = zmin + real(k - 1) * dz
            f(i,j,k) = sin(x*y*z)
          end do
        end do
      end do

      ! sample value
      x0=5.345
      y0=5.987
      z0=5.555

      interpval = interpd(Nx, Ny, Nz, dx, dy, dz, xmin, xmax,
     &                    ymin, ymax, zmin, zmax, x0, y0, z0, f)

      write(*,*) interpval

      return
      end
