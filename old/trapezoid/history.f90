subroutine history(uin,x,time,dt,dx,dy,dz,mype,npes)
   use hydro_parameters
   implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

   ! ==========================================================================
   ! Declare variables --------------------------------------------------------

   ! Arguments
   integer :: mype,npes
   real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uin
   real(dp),dimension(iu1:iu2) ::x
   real(dp) ::time,dt,dx,dy,dz

   ! Locals
   character(len=255) :: filename='history.txt'
   integer :: nguard
   integer :: l, i, j, k
   integer :: ilo, ihi, jlo, jhi, klo, khi
   integer, parameter :: LUN = 17
   real(dp) :: cell_vol
   real(dp) :: mass, energy

   ! ==========================================================================
   ! Compute the next line of the history file --------------------------------

   if (verbose) write (*,*) 'Entering history...'

   ! Some dimensional juggling ------------------------------------------------
   ! - We only want to sum over body cells, not guard cells (or we will
   !   double-count certain cells and also include the semi-artificial BCs)
   ! - The volume depends on the dimensionality
   nguard = 3

   ilo = iu1 + nguard
   ihi = iu2 - nguard
   cell_vol = dx

#if NDIM > 1
   jlo = ju1 + nguard
   jhi = ju2 - nguard
   cell_vol = cell_vol * dy
#else
   jlo = ju1
   jhi = ju2
#endif

#if NDIM > 2
   klo = ku1 + nguard
   khi = ku2 - nguard
   cell_vol = cell_vol * dz
#else
   klo = ku1
   khi = ku2
#endif

   ! Sum over local cells -----------------------------------------------------
   mass   = 0.0d0
   energy = 0.0d0
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            mass   = mass   + cell_vol * uin(1,i,j,k,1)
            energy = energy + cell_vol * uin(1,i,j,k,5)
         end do
      end do
   end do

   ! Sum over processors ------------------------------------------------------
   ! - Might we speed this up by packing all the values into an array, calling
   !   a single sumBcast, then unpacking the result?  For a larger number of
   !   variables, this may be a speed-up, but it would have to be tested.  For
   !   a small number of processors and small number of variables to
   !   communicate, it won't be significant.
   call sumBcast(mass,1)
   call sumBcast(energy,1)

   ! Print the next line of the history file ----------------------------------
   if (mype == 0) then
      if (time == 0.d0) then
         open(unit=LUN,file=filename,status='unknown',position='append')
         write(LUN,"('#',x,a13,50(3x,a13))") "time", "mass", "energy"
      else
         open(unit=LUN,file=filename,status='old',position='append')
      endif
      write(LUN,"(x,x,es13.6,50(3x,es13.6))") time, mass, energy
      close(LUN)
   endif

end subroutine history
