!----------------------------------------------------------------------
!  Debug routine
!----------------------------------------------------------------------
      real*8 function plotter(cmask,gfunc_n,g1_Np,g1_Nz,p,z,L)
        implicit none

        include 'globals.inc'
        include 'cmask.inc' 

        integer unum, stat      ! used to read/open/close a file

        integer i, j
        integer g1_Np,g1_Nz
        real*8  cmask(g1_Np,g1_Nz)
        real*8  gfunc_n(1:g1_Np,1:g1_Nz)
        real*8  p(*), z(*)
        integer L
        character(len=1024) maskfname, gfuncfname

        write (maskfname, "(A4,I1,A4)") "mask", L, ".dat"
        maskfname = trim(maskfname)
        
        write (gfuncfname, "(A5,I1,A4)") "gfunc", L, ".dat"
        gfuncfname = trim(gfuncfname)

        unum = 124
        open(unit=unum, file=maskfname, status='unknown',
     &       action='write',position='append',iostat=stat)
        open(unit=unum+1, file=gfuncfname, status='unknown',
     &       action='write',position='append',iostat=stat)
          
        do i=1,g1_Np
        do j=1,g1_Nz

          write(unum,  *, iostat=stat) p(i), z(j), cmask(i,j)

          if ( cmask(i,j) .eq. CMASK_on ) then
           write (unum+1,*, iostat=stat) p(i), z(j), gfunc_n(i,j)
          else
           write (unum+1,*, iostat=stat) p(i), z(j), -0.5d0
          endif

        end do
           write(unum,  *,iostat=stat) 
           write(unum+1,*,iostat=stat) 
        end do

        close(unit=unum,   status='keep', iostat=stat) 
        close(unit=unum+1, status='keep', iostat=stat) 

        plotter = 0d0
        return
      end
