c============================================================
c     Code to manually update some grid functions. Note
c     that RNPL supplies the routine header, declaration of
c     some routine arguments, and RETURN and END statements.
c============================================================

      integer Nx, i

      Nx = g1_Nx ! rename for convenience

      ! apply the limit at x=0
      exactK_np1(1)=1d0
      exactH_np1(1)=0d0

      do i=2,Nx
        exactK_np1(i)=x(i)/sinh(x(i))
        exactH_np1(i)=x(i)/tanh(x(i)) - 1d0
      enddo

      do i=1,Nx
        errK_np1(i)=K_np1(i)-exactK_np1(i)
        errH_np1(i)=H_np1(i)-exactH_np1(i)
      enddo

