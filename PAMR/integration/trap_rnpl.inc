c============================================================
c     Code to manually update (integrate) the grid function.
c     Note that RNPL supplies the routine header, declaration
c     of some routine arguments, and RETURN / END statements.
c============================================================

      integer i, j
      real*8  mysum     

      mysum = 0.0d0

      do i=1, g1_Np-1
      do j=1, g1_Nz-1
        mysum=mysum+dp*dz/4.0d0*
c     &                  (  p(i)   * gfunc_n(i,j)
c     &                   + p(i+1) * gfunc_n(i+1,j)
c     &                   + p(i)   * gfunc_n(i,j+1)
c     &                   + p(i+1) * gfunc_n(i+1,j+1) )
     &                 (            gfunc_n(i,j)
     &                  +           gfunc_n(i+1,j)
     &                  +           gfunc_n(i,j+1)
     &                  +           gfunc_n(i+1,j+1) )
      end do
      end do

      gfunc_np1(:,:)=gfunc_n(:,:)
c      write(*,*) 2d0*3.14159265358979d0*mysum
