cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     interp: performs bicubic interpolation, borrowing heavily from
c             Press et al, 'Numerical Recipes in Fortran 77', Ch. 3.6
c             (Cambridge University Press, 1992)
c
c             Input quantities are:
c             y   - array of length 4 containing the function values at
c                   the corner grid points of the rectangular cell in
c                   which the interpolated value is requested (counting
c                   ccw from lower left grid point)
c             y1  - same as y, but for derivative values w.r.t dir 1
c             y2  - same as y, but for derivative values w.r.t dir 2
c             y12 - same as y, but for mixed derivative
c             x1l - lower coordinates of the corner point in direction 1
c             x1u - upper coordinates of the corner point in direction 1
c             x2l - lower coordinates of the corner point in direction 2
c             x2u - upper coordinates of the corner point in direction 2
c             x1, x2 - coordinates of desired point for interpolation
c
c             The interpolated function value is returned as ansy,
c             and the interpolated derivative values as ansy1 and ansy2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,
     &ansy2)
      implicit none

      real*8, dimension(4), intent(in) :: y,y1,y2,y12
      real*8, intent(in)  :: x1l,x1u,x2l,x2u,x1,x2
      real*8, intent(out) :: ansy,ansy1,ansy2
      integer  :: i
      real*8   :: t,u
      real*8, dimension(4,4) :: c

      call coeff(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)

      t=(x1-x1l)/(x1u-x1l)
      u=(x2-x2l)/(x2u-x2l)
      
      ansy=0.0
      ansy2=0.0
      ansy1=0.0

      do i=4,1,-1
        ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
        ansy2=t*ansy2+(3.0*c(i,4)*u+2.0*c(i,3))*u+c(i,2)
        ansy1=u*ansy1+(3.0*c(4,i)*t+2.0*c(3,i))*t+c(2,i)
      end do
      
      ansy1=ansy1/(x1u-x1l)
      ansy2=ansy2/(x2u-x2l)

      end subroutine interp


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     coeff: computes necessary coefficients for interp
c
c            Input quantities are:
c            y, y1, y2, y12 - same as in interp
c            d1, d2 - length of rectangular grid in directions 1, 2
c            c - a 4x4 array storing interpolation coefficients
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeff(y,y1,y2,y12,d1,d2,c)
      implicit none

      real*8, dimension(4),   intent(in)  :: y,y1,y2,y12
      real*8, dimension(4,4), intent(out) :: c
      real*8, intent(in) :: d1,d2
      real*8, dimension(16)    :: x
      real*8, dimension(16,16) :: wt

      data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,
     &8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,
     &2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,
     &2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,
     &-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,
     &-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,
     &-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

      x(1:4)=y
      x(5:8)=y1*d1
      x(9:12)=y2*d2
      x(13:16)=y12*d1*d2

      x=matmul(wt,x)
      c=reshape(x,(/4,4/),order=(/2,1/))

      end subroutine coeff
