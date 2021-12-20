c----------------------------------------------------------------------
c evolution and other numerical routines for wave
c 
c In all routines below, y and z are only defined if dim>1 and 2 
c respectively (i.e. trying to access these arrays otherwise
c will probably cause segmentation faults)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c phi_evo() performs 1 evolution time step of the wave equation,
c solving for the field at time np1
c
c box(phi) = 0
c
c ==>
c
c phi,t,t - phi,x,x - phi,y,y - phi,z,z = 0
c
c ==>
c
c (phi_np1(i,j,k)-2*phi_n(i,j,k)+phi_nm1(i,j,k))/dt^2=
c (phi_n(i+1,j,k)-2*phi_n(i,j,k)+phi_n(i-1,j,k))/dx^2 +
c (phi_n(i,j+1,k)-2*phi_n(i,j,k)+phi_n(i,j-1,k))/dy^2 +
c (phi_n(i,j,k+1)-2*phi_n(i,j,k)+phi_n(i,j,k-1))/dz^2
c
c ==>
c
c phi_np1(i,j,k) = 2*phi_n(i,j,k) - phi_nm1(i,j,k) + 
c dt^2*((phi_n(i+1,j,k)-2*phi_n(i,j,k)+phi_n(i-1,j,k))/dx^2 +
c       (phi_n(i,j+1,k)-2*phi_n(i,j,k)+phi_n(i,j-1,k))/dy^2 +
c       (phi_n(i,j,k+1)-2*phi_n(i,j,k)+phi_n(i,j,k-1))/dz^2)
c
c----------------------------------------------------------------------
        subroutine phi_evo(phi_np1,phi_n,phi_nm1,x,y,z,dt,
     &                     dim,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,dim
        real*8 phi_np1(Nx,Ny,Nz),phi_n(Nx,Ny,Nz),phi_nm1(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt

        integer i,j,k,js,je,ks,ke
        real*8 dx,dy,dz



        dx=x(2)-x(1)
        js=1
        je=1
        ks=1
        ke=1
        if (dim.gt.1) then
           dy=y(2)-y(1)
           js=2
           je=Ny-1
           if (dim.gt.2) then
              ks=2
              ke=Nz-1
              dz=z(2)-z(1)
           end if
        end if

        do i=2,Nx-1
         do j=js,je
          do k=ks,ke

           phi_np1(i,j,k)=2*phi_n(i,j,k)-phi_nm1(i,j,k)+dt**2*
     &        (phi_n(i+1,j,k)-2*phi_n(i,j,k)+phi_n(i-1,j,k))/dx**2
           if (dim.gt.1) then
            phi_np1(i,j,k)=phi_np1(i,j,k)+dt**2*
     &         (phi_n(i,j+1,k)-2*phi_n(i,j,k)+phi_n(i,j-1,k))/dy**2 
            if (dim.gt.2) then
             phi_np1(i,j,k)=phi_np1(i,j,k)+dt**2*
     &          (phi_n(i,j,k+1)-2*phi_n(i,j,k)+phi_n(i,j,k-1))/dz**2
            end if
           end if
 
          end do
         end do
        end do

        return
        end

c----------------------------------------------------------------------
c initializes f with a generalized gaussian 'kink':
c
c f = A * exp (- (r-r0)^2/delta^2) , r>r0
c   = A , r < r0
c
c where r = sqrt ( (1-ex^2)*(xb)^2 + (1-ey^2)*(yb)^2 + (1-ez^2)*(zb)^2 )
c
c and xb = x-xu0, yb = ...
c----------------------------------------------------------------------
        subroutine gauss3d(f,A,r0,delta,xu0,yu0,zu0,ex,ey,ez,
     &                     x,y,z,dim,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,dim
        real*8 f(Nx,Ny,Nz),x(Nx),y(Ny),z(Nz)
        real*8 A,r0,delta,ex,ey,ez,xu0,yu0,zu0

        integer i,j,k
        real*8 r,xb,yb,zb

        do i=1,Nx
           xb=x(i)-xu0
           do j=1,Ny
              if (dim.lt.2) then
                 yb=0
              else
                 yb=y(j)-yu0
              end if
              do k=1,Nz

                 if (dim.lt.3) then
                    zb=0
                 else
                    zb=z(k)-zu0
                 end if
               
                 r=sqrt((1-ex**2)*xb**2+
     &                  (1-ey**2)*yb**2+
     &                  (1-ez**2)*zb**2)
                 if (r.gt.r0) then
                    f(i,j,k)=A*exp(-((r-r0)/delta)**2) 
                 else
                    f(i,j,k)=A
                 end if

              end do
           end do
        end do

        return
        end

