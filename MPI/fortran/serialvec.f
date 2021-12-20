c     This is a serial program to sum the elements of a vector

      program serialvec
      implicit none

      integer   i,     num_rows,      max_rows
      parameter ( max_rows = 100000000 )
      real*8    vect(max_rows)  ! create a real vector with max_rows entries
      real*8    summ

      num_rows = 1000
      if ( num_rows .gt. max_rows) stop "Too many numbers."

c     Create the vector
      do i = 1, num_rows
         vect(i) = float(i) ! converts i to float/real
      end do

c     Compute the sum
      summ = 0d0
      do i = 1, num_rows
         summ = summ + vect(i)
      end do

      write(*,*) summ

      stop
      end
