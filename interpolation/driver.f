ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     driver.f -- driver program for interpolation function neville()
c                 defined in interp.f. This script reads data from a 
c                 data file, picks out 5 data points surrounding a
c                 user-specified point, and then passes the desired 
c                 value and arrays to neville()
c
c                 NOTE: this program assumes the input data is
c                       specified at uniformly-spaced intervals
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program driver

      implicit none

      real*8, allocatable :: axis(:)     
      real*8, allocatable :: inputdata(:)
      real*8, dimension(0:4) :: x = 0d0, f = 0d0
      real*8  :: r, datval, dr, p
      real*8  :: neville
      integer :: i, length, interpindex, stat

      ! get length of input data
      length = 0
      open(unit=95, file='testdata', status='old', action='read')
      do
        read(95,*,iostat=stat)
        if (stat .ne. 0) exit
        length = length + 1
      end do
      close(unit=95, status='keep')

      allocate(axis(length))
      allocate(inputdata(length))

      ! read data values into arrays
      open(unit=95, file='testdata', status='old', action='read')
      do i=1,length
        read (95,*,iostat=stat) r, datval
        axis(i)=r
        inputdata(i)=datval
      end do
      close(unit=95, status='keep')
 
      ! pick the point at which you would like an interpolated value
      p=24.0002d0

      ! get array index closest to p (ASSUMES UNIFORMLY-SPACED DATA)
      dr=axis(2)-axis(1)
      interpindex=idnint(p/dr)+1
      if ( interpindex+2 .ge. length ) then
        STOP 'Interpolation error (length exceeded)'
      end if
      if ( interpindex .le. 2 ) then  ! fix for interpolation near lower bound
        interpindex = 3
      end if

      ! set appropriate array values and pass to interpolation routine
      x(0)=axis(interpindex-2)
      x(1)=axis(interpindex-1) 
      x(2)=axis(interpindex)
      x(3)=axis(interpindex+1)
      x(4)=axis(interpindex+2)

      f(0)=inputdata(interpindex-2)
      f(1)=inputdata(interpindex-1) 
      f(2)=inputdata(interpindex)
      f(3)=inputdata(interpindex+1)
      f(4)=inputdata(interpindex+2)

      write(*,*) neville(x,f,p)

      return
      end
