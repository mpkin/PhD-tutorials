!----------------------------------------------------------------------
!  A Fortran function designed to be called by a PAMR hook function.
!  Integrates gridfunction 'gfunc' over the entire domain using
!  an O(h^2) trapezoidal rule.
!----------------------------------------------------------------------
      real*8 function gfuncintpamr(cmask,gfunc_np1,gfunc_n,
     &                             g1_Np,g1_Nz,p,z,dp,dz,gw)
        implicit none

        include 'globals.inc'
        include 'cmask.inc'

        integer i, j
        integer g1_Np,g1_Nz
        real*8  cmask(g1_Np,g1_Nz)
        real*8  gfunc_np1(g1_Np,g1_Nz)
        real*8  gfunc_n(g1_Np,g1_Nz)
        real*8  dp,dz
        real*8  p(*), z(*)
        integer gw(4)
        real*8  mysum

        mysum = 0.0d0
        
        write(*,*)

        do i=1+gw(1), g1_Np-1-gw(2)
        do j=1+gw(3), g1_Nz-1-gw(4)

          if ( ( cmask(i,j)     .eq. CMASK_OFF ) .and.
     &         ( cmask(i,j+1)   .ne. CMASK_ON )  .and.
     &         ( cmask(i+1,j)   .ne. CMASK_ON )  .and.
     &         ( cmask(i+1,j+1) .ne. CMASK_ON ) ) then
          mysum=mysum+dp*dz/4.0d0*
c     &              (  p(i)   * gfunc_n(i,j)
c     &               + p(i+1) * gfunc_n(i+1,j)
c     &               + p(i)   * gfunc_n(i,j+1)
c     &               + p(i+1) * gfunc_n(i+1,j+1) )
     &              (           gfunc_n(i,j)
     &               +          gfunc_n(i+1,j)
     &               +          gfunc_n(i,j+1)
     &               +          gfunc_n(i+1,j+1) )
          end if
          
c        write(*,*) p(i), z(j), cmask(i,j)

        end do
        end do

c        gfuncintpamr=2d0*3.14159265358979d0*mysum
        gfuncintpamr=mysum

        return
      end

