      module triangle_mod
      implicit none
      private

      public :: hypo

      contains

        real function hypo(side_1,side_2)

          implicit none
 
          ! Declare calling parameters:
          real, intent(in) :: side_1      ! Length of side 1
          real, intent(in) :: side_2      ! Length of side 2

          ! Declare local variables
          real :: temp                    ! Temporary storage variable

          ! Calculate hypotenuse
          temp = side_1**2 + side_2**2
          hypo = sqrt(temp)

        end function hypo

      end module triangle_mod
