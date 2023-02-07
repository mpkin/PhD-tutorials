c-----------------------------------------------------------------------
      module globals
      implicit none

      integer counter
      real*8 Q, sigma, pi
      parameter ( Q = 1.0d0,   sigma = 1.0d0 )
      parameter ( pi = 3.14159265358979d0 )

      end module
c----------------------------------------------------------------------


c-----------------------------------------------------------------------
c
c     initdata: initializer for elliptic constraint solver
c
c-----------------------------------------------------------------------
      subroutine initdata(V,p,z,Np,Nz)
      use globals
      implicit none

      integer Np,Nz
      real*8 V(Np,Nz), p(Np), z(Nz)
      real*8 rand

      integer i, j
      real*8 Vx(Np,Nz)
      character (len=16) :: id

      counter = 0

      ! Initial guess
      do i=1,Np
      do j=1,Nz
       call RANDOM_NUMBER(rand)
       V(i,j)=rand/4d0
      end do
      end do

      return
      end subroutine


c-----------------------------------------------------------------------
c
c     lop: differential operator L(V)=0 computed where cmask=CMASK_ON
c
c-----------------------------------------------------------------------
      subroutine lop(LV,V,cmask,p,z,Np,Nz)
      use globals
      implicit none

      integer Np, Nz
      real*8 V(Np,Nz)
      real*8 cmask(Np,Nz), LV(Np,Nz)
      real*8 p(Np), z(Nz)

      integer i,j
      real*8 f(Np,Nz)
      real*8 dp, dz

      include 'cmask.inc'

      dp=(p(2)-p(1))
      dz=(z(2)-z(1))

      ! Source term
      do i=1,Np
      do j=1,Nz
        f(i,j)=Q/(sigma*sqrt(2d0*pi))**3d0
     &  *exp(-sqrt(p(i)**2d0+z(j)**2d0)**2d0
     &       /(2d0*sigma**2d0))
      end do
      end do

      do i=2,Np-1
       do j=2,Nz-1
        if (cmask(i,j).eq.CMASK_ON) then
         LV(i,j)=
     &0.1D1 / p(i) * (V(i + 1,j) - V(i - 1,j)) / dp / 0.2D1 + (V(
     &i + 1,j) - 0.2D1 * V(i,j) + V(i - 1,j)) / dp ** 2 + (V(i,j + 1) -
     &0.2D1 * V(i,j) + V(i,j - 1)) / dz ** 2 - f(i,j)
        end if
       end do
      end do

      write(*,*) counter, "lop:     ", Np, Nz

      return
      end subroutine


c-----------------------------------------------------------------------
c
c     residual: residual L[V]-rhs computed where cmask=CMASK_ON 
c
c-----------------------------------------------------------------------
      subroutine residual(res,rhs,V,cmask,p,z,norm,Np,Nz)
      use globals
      implicit none

      integer Np, Nz
      real*8 V(Np,Nz)
      real*8 cmask(Np,Nz), res(Np,Nz), rhs(Np,Nz)
      real*8 p(Np), z(Nz), norm

      integer i, j, sum

      include 'cmask.inc'

      call lop(res,V,cmask,p,z,Np,Nz)

      norm = 0d0
      sum  = 0

      do i=2,Np-1
       do j=2,Nz-1
        if (cmask(i,j).eq.CMASK_ON) then
         res(i,j) = res(i,j)-rhs(i,j)
         norm     = norm+res(i,j)**2
         sum      = sum+1
        end if
       end do
      end do

      norm = sqrt(norm/sum)

      write(*,*) counter, "residual:", Np, Nz, norm

      return
      end subroutine


c-----------------------------------------------------------------------
c
c     relax: applies alternating-direction zebra line relaxation using
c            Lapack 'dgtsv' solver
c
c-----------------------------------------------------------------------
      subroutine relax(V,V_rhs,cmask,phys_bdy,p,z,norm,Np,Nz)
      use globals
      implicit none

      integer Np, Nz
      real*8 V(Np,Nz)
      real*8 cmask(Np,Nz), V_rhs(Np,Nz)
      real*8 p(Np), z(Nz), norm
      integer phys_bdy(4)

      real*8 f(Np,Nz), Vx(Np,Nz), Vd(Np,Nz)
      real*8 dp, dz
      integer i, j, k, l
      character (len=16) :: id = 'approx'

      real*8, dimension (Nz) ::  d, du, dl, rhs
      integer nrhs, info
      real*8 res
      integer sum

      include 'cmask.inc'

      dp = (p(2)-p(1))
      dz = (z(2)-z(1))

      norm = 0d0
      sum  = 0

      ! Source term
      do i=1,Np
      do j=1,Nz
        f(i,j) = Q/(sigma*sqrt(2d0*pi))**3d0
     &  *exp(-sqrt(p(i)**2d0+z(j)**2d0)**2d0
     &       /(2d0*sigma**2d0))
      end do
      end do

      ! Exact solution
      do i=1,Np
      do j=1,Nz
        if ((i .EQ. 1) .AND. (j .EQ. Np)) then
          Vx(i,j) = -sqrt(0.2D1)*Q*pi**(-0.3D1/0.2D1)/sigma/0.4D1
        else
          Vx(i,j) = -1d0/(4d0*pi)*Q/sqrt(p(i)**2d0+z(j)**2d0)
     &    *erf(sqrt(p(i)**2d0+z(j)**2d0)/(sqrt(2d0)*sigma))
        end if
      end do
      end do

      ! Boundary condition @ p=0
      do j=1,Nz
        V(1,j) = Vx(1,j)
      end do

      ! Boundary condition @ p=infty
      do j=1,Nz
        V(Np,j) = Vx(Np,j)
      end do

      ! Boundary condition @ z=infty
      do i=1,Np
        V(i,Nz) = Vx(i,Nz)
      end do

      ! Boundary condition @ z=-infty
      do i=1,Np
        V(i,1) = Vx(i,1)
      end do

      ! Compute difference from exact solution
      do i=1,Np
      do j=1,Nz
       Vd(i,j)=V(i,j)-Vx(i,j)
      end do
      end do

      ! Output values for the initial time
      if (counter .EQ. 0) then
        id='exact'
        call w2f(Vx,Np,Nz,p,z,id)
        id='diff'
        call w2f(Vd,Np,Nz,p,z,id)
        id='approx'
        call w2f(V,Np,Nz,p,z,id)
        counter = 1
      end if

      ! X-DIRECTION
      ! Set up tridiagonal system. Note that indexing on 
      ! lower diagonal is always (j-1) when implementing the 
      ! jth equation
      do l=0,1  ! zebra sweep
      do i=2+l,Np-1-l,2
        d(1)=1.0d0
        du(1)=0.0d0
        rhs(1)=V(i,1)
        do j=2,Nz-1   ! using 5-point stencil:
          dl(j-1) =  0.1D1 / dz ** 2
          d(j)    = -2d0 / dp ** 2d0 - 2d0 / dz ** 2d0
          du(j)   = 0.1D1 / dz ** 2
          rhs(j)  = f(i,j) + V_rhs(i,j)
     #              - (V(i + 1,j)-V(i - 1,j))/p(i)/dp/0.2D1-(V
     #              (i + 1,j) + V(i - 1,j)) / dp ** 2
        end do
        dl(Nz-1) = 0.0d0
        d(Nz)    = 1.0d0
        rhs(Nz)  = V(i,Nz)
  
        nrhs = 1
        call dgtsv( Nz, nrhs, dl, d, du, rhs, Nz, info )
  
        do j=1,Nz
          V(i,j)=rhs(j)
        end do
  
      end do
      end do

      ! Y-DIRECTION
      ! Set up tridiagonal system. Note that indexing on 
      ! lower diagonal is always (j-1) when implementing the 
      ! jth equation
      do l=0,1  ! zebra sweep
      do j=2+l,Nz-1-l,2
        d(1)=1.0d0
        du(1)=0.0d0
        rhs(1)=V(1,j)
        do i=2,Np-1   ! using 5-point stencil:
          dl(i-1) = -0.1D1 / p(i) / dp / 0.2D1 + 0.1D1 / dp ** 2
          d(i)    = -2d0 / dp ** 2d0 - 2d0 / dz ** 2d0
          du(i)   = 0.1D1 / p(i) / dp / 0.2D1 + 0.1D1 / dp ** 2
          rhs(i)  =  -(V(int(i),j + 1) + V(int(i),j - 1))
     #               / dz ** 2 + f(int(i),j) + V_rhs(i,j)
        end do
        dl(Np-1) = 0.0d0
        d(Np)    = 1.0d0
        rhs(Np)  = V(Np,j)
  
        nrhs = 1
        call dgtsv( Np, nrhs, dl, d, du, rhs, Np, info )
  
        do i=1,Np
          V(i,j)=rhs(i)
        end do
  
      end do
      end do

      ! Compute residual
      do i=2,Np-1
      do j=2,Nz-1
       if (cmask(i,j).eq.CMASK_ON) then
        res =
     &0.1D1 / p(i) * (V(i + 1,j) - V(i - 1,j)) / dp / 0.2D1 + (V(
     &i + 1,j) - 0.2D1 * V(i,j) + V(i - 1,j)) / dp ** 2 + (V(i,j + 1) -
     &0.2D1 * V(i,j) + V(i,j - 1)) / dz ** 2 - f(i,j) - V_rhs(i,j)
        norm = norm + res**2
        sum  = sum + 1
       end if
      end do
      end do
      
      norm = sqrt(norm/sum)

      ! Compute difference from exact solution
      do i=1,Np
      do j=1,Nz
       Vd(i,j)=V(i,j)-Vx(i,j)
      end do
      end do

      id='exact'
      call w2f(Vx,Np,Nz,p,z,id)
      id='diff'
      call w2f(Vd,Np,Nz,p,z,id)
      id='approx'
      call w2f(V,Np,Nz,p,z,id)
      counter = counter + 1

      write(*,*) counter, "relax:   ", Np, Nz, norm

      return
      end subroutine


c-----------------------------------------------------------------------
c
c     w2f: write relevant data to text file
c
c-----------------------------------------------------------------------
      subroutine w2f(u,Np,Nz,p,z,id)
      use globals
      implicit none

      integer, intent(in) :: Np, Nz
      real*8, dimension (Np, Nz), intent(in) ::  u
      real*8, dimension (Np), intent(in) :: p
      real*8, dimension (Nz), intent(in) :: z
      character (len=16), intent(in) :: id

      character (len=16) :: num
      character (len=32) :: filename1
      integer unum, k, l, m

      write(num,'(I0)') counter
      filename1='data/sol_'//trim(id)//num

      unum = 20
 400  format( F15.10, "  ", F15.10, "  ", F16.8 )
      open(unit=unum,   file=filename1, status='replace',
     &     action='write')
      do k=1,Np
      do m=1,Nz
        write(unum,400)   p(k), z(m), u(k,m)
      end do
      end do
      close(unit=unum, status='keep')

      return
      end subroutine

