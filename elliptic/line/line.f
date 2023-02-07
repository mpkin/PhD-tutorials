c=======================================================================
c
c     line: solves the Poisson equation using alternating zebra
c           Gauss-Seidel relaxation in cylindrical coordinates
c           with fixed (exact) boundary conditions. The problem
c           solved is the Poisson equation (assuming axisymmetry):
c
c                                              2 
c             1   d   /    d          \       d
c            --- --- |  p ---  u(p,z)  |  +  --- u(p,z)  =  f(p,z)
c             p  dp   \   dp          /        2
c                                            dz 
c
c           where the source function f(p,z) is
c
c
c                              Q               / -1  p^2+z^2 \
c            f(p,z) = -------------------- exp|  --- -------  |
c                     sigma^3 sqrt(2*pi)^3     \  2  sigma^2 /
c
c
c                which allows the exact solution
c
c                       1     Q        / sqrt(p^2+z^2) \
c            u(p,z) = ---- ------- erf| --------------- |
c                     4*pi p^2+z^2     \ sigma*sqrt(2) /
c
c
c=======================================================================

      program line

      implicit none

      real*8          pmin, pmax, zmin, zmax
      parameter     ( pmin = 0.0d0,   pmax = 2.0d0 )
      parameter     ( zmin = -2.0d0,  zmax = 2.0d0 )

c     Define maximum problem size
      integer         maxnm
      parameter     ( maxnm = 2 048 )

      real*8, dimension (maxnm)        :: p, z
      real*8, dimension (maxnm, maxnm) :: u, ux, f

      integer         level
      integer         n, m, i, j
      real*8          dp, dz

      real*8          Q, sigma, pi
      parameter     ( Q = 1.0d0,   sigma = 1.0d0 )
      parameter     ( pi = 3.14159265358979d0 )

      integer         iter

      character (len=16) :: id = ''

      level = 6
      n = 2 ** level + 1
      m = 2*(2 ** level) + 1
      if( m .gt. maxnm ) then
         write(0,*) 'Insufficient internal storage'
         stop
      end if

      dp    =  (pmax - pmin) / (n - 1)
      dz    =  (zmax - zmin) / (m - 1)

      do i = 1, n
         p(i) = pmin + (i-1) * dp
      end do

      do j = 1, m
         z(j) = zmin + (j-1) * dz
      end do

      do i=1,n
      do j=1,m
        f(i,j)=Q/(sigma*sqrt(2d0*pi))**3d0
     &  *exp(-sqrt(p(i)**2d0+z(j)**2d0)**2d0
     &       /(2d0*sigma**2d0))
      end do
      end do

      do i=1,n
      do j=1,m
        if ((i .EQ. 1) .AND. (j .EQ. n)) then
          ux(i,j)=-sqrt(0.2D1)*Q*pi**(-0.3D1/0.2D1)/sigma/0.4D1
        else
          ux(i,j)=-1d0/(4d0*pi)*Q/sqrt(p(i)**2d0+z(j)**2d0)
     &    *erf(sqrt(p(i)**2d0+z(j)**2d0)/(sqrt(2d0)*sigma))
        end if
      end do
      end do

c     Set the initial guess
      do i=1,n
      do j=1,m
        u(i,j)=0d0
      end do
      end do

c     Boundary condition @ p=0
      do j=1,m
        u(1,j) = ux(1,j)
      end do

c     Boundary condition @ p=infty
      do j=1,m
      u(n,j) = ux(n,j)
      end do

c     Boundary condition @ z=infty
      do i=1,n
      u(i,m) = ux(i,m)
      end do

c     Boundary condition @ z=-infty
      do i=1,n
      u(i,1) = ux(i,1)
      end do

      do iter=1,5000
        write(*,*) 'iter = ', iter
        call relax(u,f,maxnm,n,m,dp,dz,p,z)
      end do

      call w2f(u,ux,f,maxnm,n,m,p,z,id)

      stop
      end

c-----------------------------------------------------------------------
c
c     relax: solves a 1D BVP using Gauss-Seidel relaxation via LAPACK,
c     assuming a 5-point stencil for the Laplacian
c
c-----------------------------------------------------------------------
      subroutine relax(u,f,maxnm,n,m,dp,dz,p,z)

      implicit none

      real*8, dimension (maxnm, maxnm), intent(inout) :: u

      integer, intent(in) :: maxnm, n, m
      real*8, dimension (maxnm, maxnm), intent(in) :: f
      real*8, intent(in) :: dp, dz
      real*8, dimension (maxnm), intent(in) :: p, z

      real*8, dimension (maxnm) :: d, du, dl, rhs
      integer :: i, j, l
      integer :: nrhs, info

c     X-DIRECTION
c     Set up tridiagonal system. Note that indexing on 
c     lower diagonal is always (j-1) when implementing the 
c     jth equation
      do l=0,1  ! zebra sweep
      do i=2+l,n-1-l,2
        d(1)=1.0d0
        du(1)=0.0d0
        rhs(1)=u(i,1)
        do j=2,m-1   ! using 5-point stencil:
          dl(j-1) =  0.1D1 / dz ** 2
          d(j)    = -2d0 / dp ** 2d0 - 2d0 / dz ** 2d0
          du(j)   = 0.1D1 / dz ** 2
          rhs(j)  = f(i,j) - (u(i + 1,j)-u(i - 1,j))/p(i)/dp/0.2D1-(u
     #              (i + 1,j) + u(i - 1,j)) / dp ** 2
        end do
        dl(m-1) = 0.0d0
        d(m)    = 1.0d0
        rhs(m)  = u(i,m)
  
        nrhs = 1
        call dgtsv( m, nrhs, dl, d, du, rhs, m, info )
  
        do j=1,m
          u(i,j)=rhs(j)
        end do
  
      end do
      end do

c     Y-DIRECTION
c     Set up tridiagonal system. Note that indexing on 
c     lower diagonal is always (j-1) when implementing the 
c     jth equation
      do l=0,1  ! zebra sweep
      do j=2+l,m-1-l,2
        d(1)=1.0d0
        du(1)=0.0d0
        rhs(1)=u(1,j)
        do i=2,n-1   ! using 5-point stencil:
          dl(i-1) = -0.1D1 / p(i) / dp / 0.2D1 + 0.1D1 / dp ** 2
          d(i)    = -2d0 / dp ** 2d0 - 2d0 / dz ** 2d0
          du(i)   = 0.1D1 / p(i) / dp / 0.2D1 + 0.1D1 / dp ** 2
          rhs(i)  =  -(u(int(i),j + 1) + u(int(i),j - 1))
     #               / dz ** 2 + f(int(i),j)
        end do
        dl(n-1) = 0.0d0
        d(n)    = 1.0d0
        rhs(n)  = u(n,j)
  
        nrhs = 1
        call dgtsv( n, nrhs, dl, d, du, rhs, n, info )
  
        do i=1,n
          u(i,j)=rhs(i)
        end do
  
      end do
      end do

      return
      end subroutine


c-----------------------------------------------------------------------
c
c     w2f: write relevant data to text file
c
c-----------------------------------------------------------------------
      subroutine w2f(u,ux,f,maxnm,n,m,p,z,id)

      implicit none

      integer, intent(in) :: maxnm, n, m
      real*8, dimension (maxnm, maxnm), intent(in) :: u, ux, f
      real*8, dimension (maxnm), intent(in) :: p, z
      character (len=16), intent(in) :: id

      integer :: unum, k, l
      character (len=32) :: filename1, filename2, filename3, filename4
      
      filename1='data/sol_numerical'//id
      filename2='data/sol_exact'//id
      filename3='data/sol_diff'//id
      filename4='data/sol_source'//id

      unum = 20
 400  format( F15.10, "  ", F15.10, "  ", F12.7 )
      open(unit=unum,   file=filename1, status='replace',
     &     action='write')
      open(unit=unum+1, file=filename2, status='replace',
     &     action='write')
      open(unit=unum+2, file=filename3, status='replace',
     &     action='write')
      open(unit=unum+3, file=filename4, status='replace',
     &     action='write')
      do k=1,n
      do l=1,m
        write(unum,400)   p(k), z(l), u(k,l)
        write(unum+1,400) p(k), z(l), ux(k,l)
        write(unum+2,400) p(k), z(l), u(k,l)-ux(k,l)
        write(unum+3,400) p(k), z(l), f(k,l)
      end do
      end do
      close(unit=unum, status='keep')
      close(unit=unum+1, status='keep')
      close(unit=unum+2, status='keep')
      close(unit=unum+3, status='keep')

      return
      end subroutine
