c=======================================================================
c
c     solver: solves the Poisson equation in cylindrical coordinates
c             using successive over-relaxation (SOR)
c
c                                              2 
c             1   d   /    d          \       d
c            --- --- |  p ---  U(p,z)  |  +  --- U(p,z)  =  0
c             p  dp   \   dp          /        2
c                                            dz 
c
c=======================================================================

      program solver

      implicit none

      real*8        pi
      parameter   ( pi=3.14159265358979 )

c     2D allocatable arrays to store the grid values and source
      real*8, allocatable, dimension(:,:) ::    u,  uold, f

c     2D allocatable arrays to store cylindrical coordinate variables
      real*8, allocatable, dimension(:) ::    p,  z

      integer       i,  j,  k       ! loop variables
      real*8        x,  y           ! rectangular coordinate variables

      integer       Np, Nz          ! grid size [Np x Nz] for (p, z)
      real*8        dp, dz          ! discretization scales
      real*8        pmin, zmin      ! solution domain
      real*8        pmax, zmax      ! solution domain
      real*8        Q,  sigma       ! Gaussian potential parameters
      real*8        exp_arg         ! used for potential calculation
      real*8        w               ! overrelaxation parameter
      integer       level, stepsize ! output level/stepsize
      integer       iter            ! current iteration
      integer       maxiter         ! max # of iterations
      real*8        err             ! l-infinity norm
      real*8        tol             ! error tol (stopping criteria)

      integer       unum, stat      ! used to read/open/close a file

c     Read initial parameters from a file. The parameters must be
c     separated by a newline and must be organized in the order:
c     Np        Grid size in p
c     Nz        Grid size in z
c     pmax      Solution domain (p,z) in [0,pmax]x[-zmax,zmax]
c     zmax      Solution domain (p,z) in [0,pmax]x[-zmax,zmax]
c     sigma     Width of the source potential (Gaussian)
c     Q         Charge (i.e. height) of the source potential (Gaussian)
c     tol       Error tol (stopping criteria) for iteration algorithm
c     maxiter   Maximum allowed iterations before forced exit
c     w         Overrelaxation parameter
c     level     Output level (for convergence testing)

      unum = 114
      open(unit=unum, file='initial_params', status='old',
     &     action='read',iostat=stat)
      read(unum,*) Np
      read(unum,*) Nz
      read(unum,*) pmax
      read(unum,*) zmax 
      read(unum,*) sigma
      read(unum,*) Q
      read(unum,*) tol
      read(unum,*) maxiter
      read(unum,*) w
      read(unum,*) level
      close(unit=unum, status='keep', iostat=stat)

      pmin=0d0
      zmin=-1d0*zmax

      dp=abs((pmax-pmin)/real(Np))
      dz=abs((zmax-zmin)/real(Nz))
      
      write(*,*) dp, dz

c     Allocate the solution arrays
      allocate( u(1:Np,1:Nz) )
      allocate( uold(1:Np,1:Nz) )  
      allocate( f(1:Np,1:Nz) )
      
c     Fill the coordinate arrays
      allocate( p(1:Np) )
      allocate( z(1:Nz) )
      p(1)=0d0
      z(1)=-1d0*zmax
      do i=2,Np
        p(i)=p(i-1)+dp
      enddo
      do j=2,Nz
        z(j)=z(j-1)+dz
      enddo

      write(*,*)
      write(*,101) maxiter
      write(*,102) tol 
      write(*,103) w 
      write(*,*)
 101  format( 'Maximum iterations:        ', I13 )
 102  format( 'Specified tolerance:       ', E13.8 )
 103  format( 'Overrelaxation parameter:  ', E13.8 )

c     Set the initial guess
      do i=1,Np
      do j=1,Nz
        u(i,j)=1.0  ! initial guess is constant everywhere
        f(i,j)=0.0  
        exp_arg=-1.0*(p(i)**2.0+z(j)**2.0)/(2.0*(sigma**2.0))
        if ( exp_arg .LT. -500.0 ) then
          f(i,j)=0.0  ! bug fix for exp() underflow
        else
          f(i,j)=(Q/(sigma**3.0*sqrt(2*pi)**3.0))
     &            *exp(exp_arg) 
        endif
      end do
      end do


c     Begin main solver loop
      do iter=1,maxiter
      uold=u
      err = 0d0   

c     boundary condition at p=0 (i.e. i=1)
      i=1
      do j=2,Nz-1

      ! Quadratic fit
      u(i,j) = 0.4D1 / 0.3D1 * u(2,int(j)) - u(3,int(j)) / 0.3D1

      end do

c     boundary condition at p=pmax (i.e. i=Np)
      i=Np
      do j=1,Nz

      u(i,j)=-(-4* u(Np - 1,j) + u(Np - 2,j)) / (2*dp*p(i) + 3*p(i)
     # ** 2 + 3 * z(j) ** 2) * (p(i) ** 2 + z(j) ** 2)

      end do
      
c     interior points (copied over from previous versions of the code)
      do i=2,Np-1
      do j=2,Nz-1

      u(i,j)=u(i,j) - w*(2 * dp ** 2 * dz ** 2*p(i)* f(i,j) - 2 * dp **
     & 2*p(i)* u(i,j - 1) - 2 * dp ** 2*p(i)* u(i,j + 1) + dp * dz ** 2 
     &* u(i - 1,j) - dp * dz ** 2 * u(i + 1,j) - 2 * dz ** 2*p(i)* u(i -
     & 1,j) - 2 * dz ** 2 *p(i)*u(i + 1,j) + 4 * u(i,j)*p(i)* (dp ** 2 +
     & dz ** 2)) /p(i)/(dp ** 2 + dz ** 2)

      end do 
      end do

c     boundary condition at z=zmin (i.e. j=1)
      j=1
      do i=1,Np

      u(i,j)=-(4*u(i,2) - u(i,3)) / (2 * dz * z(j) - 3*p(i)** 2 - 3 
     #*z(j)** 2) * (p(i)** 2 +z(j)** 2)

      end do

c     boundary condition at z=zmax (i.e. j=Nz)
      j=Nz
      do i=1,Np

      u(i,j)=-(-4*u(i,Nz - 1) + u(i,Nz - 2)) / (2*dz*z(j)+3*p(i)
     # ** 2 + 3 * z(j)** 2) * (p(i)** 2 + z(j)** 2)

      end do
      
c     add up the errors (l2-norm)
      do i=1,Np
      do j=1,Nz
       err = err + (uold(i,j)-u(i,j))**2.0d0
      end do
      end do

      if ( sqrt(err)/(Np*Nz) < tol ) exit
      if ( mod(iter,1000) .EQ. 0) write(*,*) iter

      end do  ! end main iterative do-loop

      if ( iter < maxiter ) then
 200   format( 'Tolerance satisfied in ', I6, ' iterations.' )
       write(*,200) iter
      else
       write(*,300) err
 300   format( 'Maxiter reached. Program exited with final error',
     &           F15.12 )
      endif

      write(*,*)

!     Compute the stepsize for outputting the data
      if ( level .eq. 0 ) then
        stepsize = 1
      elseif ( level .eq. 1 ) then
        stepsize = 2
      elseif ( level .eq. 2 ) then
        stepsize = 4
      elseif ( level .eq. 3 ) then
        stepsize = 8
      endif

!     Write to file for the solution & source in x-y coords
      unum = 110
      open(unit=unum, file='solution_data', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum+1, file='source_data', status='replace',
     &     action='write',iostat=stat)
      do i=1,Np,stepsize
       do j=1,Nz,stepsize
        write(unum,400,iostat=stat)   p(i), z(j), u(i,j)
        write(unum+1,400,iostat=stat) p(i), z(j), f(i,j)
       end do
       write(unum,*)   ! (gnuplot compatibility)
       write(unum+1,*) ! (gnuplot compatibility)
      end do
 400  format( F10.3, "  ", F10.3, "  ", F19.14 )
      close(unit=unum,   status='keep', iostat=stat)
      close(unit=unum+1, status='keep', iostat=stat)

      call sol_checker(pi,u,Np,Nz,dp,dz,pmax,zmax,Q,sigma,stepsize)

      stop
      end


c-----------------------------------------------------------------------
c
c     subroutine sol_checker: Assuming a specific potential, this
c     routine computes and outputs the exact solution and absolute error 
c     by comparing it to the numerical solution computed in the main
c     program. The potential is assumed to be
c
c                              Q               / -1  p^2+z^2 \
c            f(p,z) = -------------------- exp|  --- -------  |
c                     sigma^3 sqrt(2*pi)^3     \  2  sigma^2 /
c           
c     which has the exact solution
c           
c                       1     Q        / sqrt(p^2+z^2) \
c            u(p,z) = ---- ------- erf| --------------- |
c                     4*pi p^2+z^2     \ sigma*sqrt(2) /
c           
c     which can be confirmed by direct substitution. This example is
c     taken from
c
c!  https://en.wikipedia.org/wiki/Poisson's_equation#Potential_of_a_Gaussian_charge_density
c
c-----------------------------------------------------------------------

      subroutine sol_checker(pi,u,Np,Nz,dp,dz,pmax,zmax,Q,sigma,
     #                       stepsize)

      implicit none

      real*8        pi

      real*8, dimension(1:Np,1:Nz) ::   u
      real*8, dimension(1:Np) ::   p
      real*8, dimension(1:Nz) ::   z

      integer       Np, Nz, stepsize
      real*8        dp, dz, pmax, zmax, Q,  sigma
      
      real*8        exact,  abserr
      integer       i,  j
      integer       unum, stat    ! used to read/open/close a file

c     Fill the coordinate arrays
      p(1)=0d0
      z(1)=-1d0*zmax
      do i=2,Np
        p(i)=p(i-1)+dp
      enddo
      do j=2,Nz
        z(j)=z(j-1)+dz
      enddo


!c    Write the file for the exact solution & abserr
      unum = 112
      open(unit=unum, file='exactsol_data', status='replace',
     &     action='write',iostat=stat)
      open(unit=unum+1, file='abserr_data', status='replace',
     &     action='write',iostat=stat)
      do i=1,Np,stepsize
      do j=1,Nz,stepsize
        if (( p(i) .EQ. 0d0 ) .AND. ( z(j) .EQ. 0d0 )) then
          exact=-0.634936358   ! determined graphically (see Maple file)
        else
        exact=-(1.0/(4.0*pi))*(Q/sqrt(p(i)**2.0+z(j)**2.0))
     &        *erf(sqrt(p(i)**2.0+z(j)**2.0)/(sigma*sqrt(2.0)))
        endif
        abserr=abs(u(i,j)-exact)
        write(unum,400,iostat=stat)   p(i), z(j), exact
        write(unum+1,400,iostat=stat) p(i), z(j), abserr
      end do
      write(unum,*)   ! (gnuplot compatibility)
      write(unum+1,*) ! (gnuplot compatibility)
      end do
 400  format( F10.3, "  ", F10.3, "  ", F13.8 )
      close(unit=unum, status='keep', iostat=stat)
      close(unit=unum+1, status='keep', iostat=stat)

      return
      end
