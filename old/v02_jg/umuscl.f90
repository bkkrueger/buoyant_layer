!> \file
!! contains subroutines umuscl(), cmpflxm(), cmp_mag_flx(),
!! fast_mhd_speed() and uslope()
!<
!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine umuscl
!
!> Performs a single timestep evolution of mhd.
!<
subroutine umuscl(uin,gravin,flux,flux_pre_tot,emfx,emfy,emfz,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const            
  use hydro_parameters
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2       ,1:ndim) :: flux_pre_tot
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2              )::emfx,emfy,emfz

  ! Primitive variables
  real(dp),dimension(:,:,:,:,:),allocatable::qin 
  real(dp),dimension(:,:,:,:,:),allocatable::bf  

  ! Slopes
  real(dp),dimension(:,:,:,:,:,:),allocatable::dq

  ! Slopes of the staggered magnetic field 
  real(dp),dimension(:,:,:,:,:,:),allocatable::dbf

  ! Left and right state arrays, all values are defined on the centre of the faces
  real(dp),dimension(:,:,:,:,:,:),allocatable::qm,qp
  real(dp),dimension(:,:,:,:,:,:),allocatable::qRT,qRB,qLT,qLB

  ! Intermediate fluxes
  real(dp),dimension(:,:,:,:,:),allocatable::fx
  real(dp),dimension(:,:,:,:  ), allocatable :: fx_pre_tot

  ! emf's
  real(dp), dimension(:,:,:,:    ), allocatable :: emf

  ! Local scalar variables
  integer::i,j,k,l,ivar,idim
  integer::ilo,ihi,jlo,jhi,klo,khi

  if (verbose) write (*,*) 'Entering umuscl...'

  allocate(qin  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar       ))
  allocate(bf   (1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3    ))
  allocate(dq   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(dbf  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3   ,1:ndim))
  allocate(qm   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qp   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim))
  allocate(qRT  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3   ))
  allocate(qRB  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3   ))
  allocate(qLT  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3   ))
  allocate(qLB  (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3   ))
  allocate(fx   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar       ))

  allocate(fx_pre_tot(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2        ))

  if (ndim>1) allocate(emf (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2))

  ilo=min(1,iu1+3); ihi=max(1,iu2-3)
  jlo=min(1,ju1+3); jhi=max(1,ju2-3)
  klo=min(1,ku1+3); khi=max(1,ku2-3)

  ! Compute primitive variables
  if(verbose) write(*,*) 'Compute primitive variables...'
  call ctoprim(uin,qin,bf,gravin,dt,ngrid)

  ! Compute TVD slopes
  if(verbose) write(*,*) 'Compute TVD slopes...'
  call uslope(bf,qin,dq,dbf,ngrid)

  ! Compute 3D traced-states in all three directions
  if(verbose) write(*,*) 'Compute predicted states...'
  if (ndim==3) then
     call trace3d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid)
  else if (ndim==2) then
     call trace2d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy   ,dt,ngrid)
  else if (ndim==1) then
     call trace1d(qin,   dq,    qm,qp,                dx      ,dt,ngrid)
  endif

  !Gravity predictor step
  do idim=1,ndim
     call gravity_predictor(qm (:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
     call gravity_predictor(qp (:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
     call gravity_predictor(qRT(:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
     call gravity_predictor(qRB(:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
     call gravity_predictor(qLT(:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
     call gravity_predictor(qLB(:,:,:,:,2:4,idim),gravin,dt,.true.,.true.,.true.,ngrid)
  end do

  ! Cooling predictor step
  do idim=1,ndim
     call cooling_predictor(qm (:,:,:,:,1,idim),qm (:,:,:,:,5,idim))
     call cooling_predictor(qp (:,:,:,:,1,idim),qp (:,:,:,:,5,idim))
  end do

  ! 1D flux in X direction
  if(verbose) write(*,*) 'Compute fluxes (x-direction)...'
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,6,7,8, &
       &       fx,fx_pre_tot,ngrid)

  ! SAVE flux in output array
  do ivar=1,nvar
     do k=kf1,kf2
        do j=jf1,jf2
           do i=if1,if2
              do l=1,ngrid
                 flux        (l,i,j,k,ivar,1)=fx        (l,i,j,k,ivar)*dt
                 flux_pre_tot(l,i,j,k,     1)=fx_pre_tot(l,i,j,k     )*dt
              end do
           end do
        end do
     end do
  end do

  !2D 
  ! 1D flux in Y direction
  if (ndim>1) then
     if(verbose) write(*,*) 'Compute fluxes (y-direction)...'
     call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
          &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
          &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,7,6,8, &
          &       fx,fx_pre_tot,ngrid)
     ! Save flux in output array
     do ivar=1,nvar
        do k=kf1,kf2
           do j=jf1,jf2
              do i=if1,if2
                 do l=1,ngrid
                    flux        (l,i,j,k,ivar,2)=fx         (l,i,j,k,ivar)*dt
                    flux_pre_tot(l,i,j,k,     2)=fx_pre_tot(l,i,j,k     )*dt
                 end do
              end do
           end do
        end do
     end do

     call cmp_mag_flx(qRT,iu1+1,iu2+1,ju1+1,ju2+1,ku1  ,ku2  , &
          &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
          &           qLT,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
          &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
          &               if1  ,if2  ,jf1  ,jf2  ,klo  ,khi  , 2,3,4,6,7,8, &
          &           emf,ngrid)

     ! Save vector in output array
     do k=klo,khi
        do j=jf1,jf2
           do i=if1,if2
              do l=1,ngrid
                 emfz(l,i,j,k)=emf(l,i,j,k)*dt
              end do
           end do
        end do
     end do
  endif

  !3D
  !  1D flux in Z direction
  if (ndim>2) then
     if(verbose) write(*,*) 'Compute fluxes (z-direction)...'
     call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
          &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
          &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,8,6,7, &
          &       fx,fx_pre_tot,ngrid)
     ! Save flux in output array
     do ivar=1,nvar
        do k=kf1,kf2
           do j=jf1,jf2
              do i=if1,if2
                 do l=1,ngrid
                    flux        (l,i,j,k,ivar,3)=fx        (l,i,j,k,ivar)*dt
                    flux_pre_tot(l,i,j,k,     3)=fx_pre_tot(l,i,j,k     )*dt
                 end do
              end do
           end do
        end do
     end do
     ! Calculate Ey
     call cmp_mag_flx(qRT,iu1+1,iu2+1,ju1,ju2,ku1+1,ku2+1, &
          &           qLT,iu1  ,iu2  ,ju1,ju2,ku1+1,ku2+1, &
          &           qRB,iu1+1,iu2+1,ju1,ju2,ku1  ,ku2  , &
          &           qLB,iu1  ,iu2  ,ju1,ju2,ku1  ,ku2  , &
          &               if1  ,if2  ,jlo,jhi,kf1  ,kf2  , 4,2,3,8,6,7, &
          &           emf,ngrid)
     ! Save vector in output array
     do k=kf1,kf2
        do j=jlo,jhi
           do i=if1,if2
              do l=1,ngrid
                 emfy(l,i,j,k)=emf(l,i,j,k)*dt
              end do
           end do
        end do
     end do
     ! Calculate Ex
     call cmp_mag_flx(qRT,iu1,iu2,ju1+1,ju2+1,ku1+1,ku2+1, &
          &           qRB,iu1,iu2,ju1+1,ju2+1,ku1  ,ku2  , &
          &           qLT,iu1,iu2,ju1  ,ju2  ,ku1+1,ku2+1, &
          &           qLB,iu1,iu2,ju1  ,ju2  ,ku1  ,ku2  , &
          &               ilo,ihi,jf1  ,jf2  ,kf1  ,kf2  , 3,4,2,7,8,6, &
          &           emf,ngrid)

     ! Save vector in output array
     do k=kf1,kf2
        do j=jf1,jf2
           do i=ilo,ihi
              do l=1,ngrid
                 emfx(l,i,j,k)=emf(l,i,j,k)*dt
              end do
           end do
        end do
     end do
  endif

  deallocate(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,fx,fx_pre_tot)
  if (ndim>1) deallocate(emf)

  if (verbose) write(*,*) 'End of umuscl...'

  return

end subroutine umuscl
!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CMPFLXM
!
!> UNSPLIT      Unsplit first (second) order Godunov integrator for
!!              polytropic magnetized gas dynamics using 
!!              the Lax-Friedrish sheme
!!              The mesh for the magnetic field is staggered 
!!              The sheme follows closely the paper by
!!              Londrillo & Del Zanna ApJ 2000, 530, 508, 
!!
!!  Inputs/Outputs
!!  uin         => (const)  input state
!! \f$ u_{\mathrm{in}} = (\rho, \rho v_{x}, \rho v_{y}, \rho v_{z}, E_{\mathrm{tot}}, B_{x}, B_{y}, B_{z} )\f$
!!       the hydro variable are centered whereas \f$B_{x}, B_{y}, B_{z}\f$ are written 
!!       on the face (staggered mesh)
!!       note that in MHD one can be in 1 or 2D and have 3 components for v and B
!!       therefore one has 2 variables for the dimension : "ndim" the spatial dimension
!!       
!!  gravin      => (const)  input gravitational acceleration
!!  iu1,iu2     => (const)  first and last index of input array,
!!  ju1,ju2     => (const)  cell centered,    
!!  ku1,ku2     => (const)  including buffer cells.
!!  flux       <=  (modify) return fluxes in the 3 coord directions
!!  if1,if2     => (const)  first and last index of output array,
!!  jf1,jf2     => (const)  edge centered,
!!  kf1,kf2     => (const)  for active cells only.
!!  dx,dy,dz    => (const)  \f$(dx,dy,dz)\f$
!!  dt          => (const)  time step
!!  ngrid       => (const)  number of sub-grids
!!  ndim        => (const)  number of dimensions
!<
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi,ln,lt1,lt2,bn,bt1,bt2, &
     &             flx,flx_pre_tot,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use variables , only : x,y,dx,dy,dz,ds,dv,cartesian,cylindrical,spherical,fargo

  implicit none

  integer :: ngrid
  integer :: ln,lt1,lt2,bn,bt1,bt2
  integer :: im1,im2,jm1,jm2,km1,km2
  integer :: ip1,ip2,jp1,jp2,kp1,kp2
  integer :: ilo,ihi,jlo,jhi,klo,khi

  real(dp), dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim) :: qm
  real(dp), dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim) :: qp 
  real(dp), dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar       ) :: flx
  real(dp), dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2              ) :: flx_pre_tot

  ! local variables
  integer :: i,j,k,l,n,idim,xdim,m

  real(dp), dimension(:,:,:,:), allocatable :: rS,uS,vS,wS,AS,BS,CS,ptotS  
  real(dp), dimension(:      ), allocatable :: qleft,qright,qgdnv
  real(dp), dimension(:      ), allocatable :: fgdnv,umean
  real(dp) :: bn_mean,zero_flux,weight_l,weight_r,S_l,S_r,shear,Ekin,Emag,Etot
  real(dp) :: ro,uo,vo,wo,bo,co,Ptoto,pressure,cotanxc,rc

  allocate(qleft(nvar),qright(nvar),qgdnv(nvar))
  allocate(fgdnv(nvar),umean(nvar)             )

  allocate(rS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(uS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(vS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(wS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(AS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(BS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(CS   (1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))
  allocate(ptotS(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2))

  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Enforce continuity for normal magnetic field
              bn_mean = half*(qm(l,i,j,k,bn,xdim)+qp(l,i,j,k,bn,xdim))

              ! Left state
              qleft (1) = qm(l,i,j,k,1  ,xdim) ! Mass density
              qleft (2) = qm(l,i,j,k,5  ,xdim) ! Pressure
              qleft (3) = qm(l,i,j,k,ln ,xdim) ! Normal velocity
              qleft (4) = bn_mean            ! Normal magnetic field
              qleft (5) = qm(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qleft (6) = qm(l,i,j,k,bt1,xdim) ! Tangential magnetic field 1
              qleft (7) = qm(l,i,j,k,lt2,xdim) ! Tangential velocity 2
              qleft (8) = qm(l,i,j,k,bt2,xdim) ! Tangential magnetic field 2

              ! Right state
              qright(1) = qp(l,i,j,k,1  ,xdim) ! Mass density
              qright(2) = qp(l,i,j,k,5  ,xdim) ! Pressure
              qright(3) = qp(l,i,j,k,ln ,xdim) ! Normal velocity
              qright(4) = bn_mean              ! Normal magnetic field
              qright(5) = qp(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qright(6) = qp(l,i,j,k,bt1,xdim) ! Tangential magnetic field 1
              qright(7) = qp(l,i,j,k,lt2,xdim) ! Tangential velocity 2
              qright(8) = qp(l,i,j,k,bt2,xdim) ! Tangential magnetic field 2

              ! Other advected quantities
              do n = 9, nvar
                 qleft (n) = qm(l,i,j,k,n,xdim)
                 qright(n) = qp(l,i,j,k,n,xdim)    
              end do

              ! Solve 1D Riemann problem
              zero_flux = one
              select case (iriemann)
              case (iroe)
                 call athena_roe    (qleft,qright,fgdnv,zero_flux)
              case (illf)
                 call lax_friedrich (qleft,qright,fgdnv,zero_flux)
              case (ihll)
                 call hll           (qleft,qright,fgdnv          )
              case (ihlld)
                 call hlld          (qleft,qright,fgdnv          ,ro,uo,vo,wo,bo,co,Ptoto)
              case (iupwind)
                 call lax_friedrich (qleft,qright,fgdnv,zero_flux)
              case (iacoustic)
                 call hydro_acoustic(qleft,qright,fgdnv          )
              end select

              ! Upwind solver in case of the shearing box (only on the hydro variables))
              if (cartesian) then
#if NDIM==1
                 if (Omega0>0.d0) then
                    shear=-1.5*Omega0*(x(i)+x(i-1))
                    fgdnv(8) = fgdnv(8) + shear*bn_mean
                 endif
#endif
#if NDIM==2
                 if (Omega0>0.d0) then
                    if (xdim==1) shear=-1.5*Omega0*(x(i)+x(i-1))
                    if (xdim==2) shear=-1.5*Omega0*x(i)
                    fgdnv(8) = fgdnv(8) + shear*bn_mean
                 endif
#endif
#if NDIM==3
                 if ((Omega0>0.d0).and.(xdim==2).and..not.(fargo)) then
                    shear=-1.5*Omega0*x(i)
                    if (shear>0.d0) then
                       Emag = half*(qleft(4)**2+qleft(6)**2+qleft(8)**2)
                       Ekin = half*(qleft(3)**2+qleft(5)**2+qleft(7)**2)
                       Etot = Ekin + Emag + qleft(2)/(gamma-1.d0)
                       fgdnv(1) = fgdnv(1) + shear*qleft(1)
                       fgdnv(2) = fgdnv(2) + shear*(Etot+Emag-bn_mean**2)
                       fgdnv(3) = fgdnv(3) + shear*qleft(1)*qleft(3)
                       fgdnv(5) = fgdnv(5) + shear*qleft(1)*qleft(5)
                       fgdnv(7) = fgdnv(7) + shear*qleft(1)*qleft(7)
                    else
                       Emag = half*(qright(4)**2+qright(6)**2+qright(8)**2)
                       Ekin = half*(qright(3)**2+qright(5)**2+qright(7)**2)
                       Etot = Ekin + Emag + qright(2)/(gamma-1.d0)
                       fgdnv(1) = fgdnv(1) + shear*qright(1)
                       fgdnv(2) = fgdnv(2) + shear*(Etot+Emag-bn_mean**2)
                       fgdnv(3) = fgdnv(3) + shear*qright(1)*qright(3)
                       fgdnv(5) = fgdnv(5) + shear*qright(1)*qright(5)
                       fgdnv(7) = fgdnv(7) + shear*qright(1)*qright(7)
                    endif
                 endif
#endif
              endif

              !Swap the variables in order to match the uin convention
              flx(l,i,j,k,1)   = fgdnv(1)*ds(i,j,k,xdim) !density
              flx(l,i,j,k,5)   = fgdnv(2)*ds(i,j,k,xdim) !total energy
              flx(l,i,j,k,ln)  = fgdnv(3)*ds(i,j,k,xdim) !normal velocity
              flx(l,i,j,k,bn)  = fgdnv(4)*ds(i,j,k,xdim) !normal magnetic field
              flx(l,i,j,k,lt1) = fgdnv(5)*ds(i,j,k,xdim) !transverse velocities 1
              flx(l,i,j,k,bt1) = fgdnv(6)*ds(i,j,k,xdim) !transverse magnetic field 1
              flx(l,i,j,k,lt2) = fgdnv(7)*ds(i,j,k,xdim) !transverse velocities 2
              flx(l,i,j,k,bt2) = fgdnv(8)*ds(i,j,k,xdim) !transverse magnetic field 2

              !Store star state
              rS(l,i,j,k)    = ro
              uS(l,i,j,k)    = uo
              vS(l,i,j,k)    = vo
              wS(l,i,j,k)    = wo
              AS(l,i,j,k)    = bn_mean
              BS(l,i,j,k)    = bo
              CS(l,i,j,k)    = co
              ptotS(l,i,j,k) = Ptoto

           end do
        end do
     end do
  end do

  ! Add geometrical terms
  if (xdim==1) then
     if (cartesian) then
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi-1
                 do l = 1, ngrid
                    flx_pre_tot(l,i,j,k)=-(ptotS(l,i+1,j,k)*ds(i+1,j,k,1)-ptotS(l,i,j,k)*ds(i,j,k,1))
                 end do
              end do
           end do
        end do
     endif
     if (cylindrical) then
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi-1              
                 do l = 1, ngrid
                    flx_pre_tot(l,i,j,k)=-dv(i,j,k)*(ptotS(l,i+1,j,k)-ptotS(l,i,j,k))/dx
                    ro=half*(rS(l,i,j,k)+rS(l,i+1,j,k))
                    vo=half*(vS(l,i,j,k)+vS(l,i+1,j,k))
                    bo=half*(BS(l,i,j,k)+BS(l,i+1,j,k))
                    !pressure=(ro*vo*vo-bo*bo)*(ds(i+1,j,k,1)-ds(i,j,k,1))
                    pressure=-bo*bo*(ds(i+1,j,k,1)-ds(i,j,k,1))
                    flx_pre_tot(l,i,j,k)=flx_pre_tot(l,i,j,k) + pressure
                 end do
              end do
           end do
        end do
     endif
     if (spherical) then
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi-1              
                 do l = 1, ngrid
                    flx_pre_tot(l,i,j,k)=-dv(i,j,k)*(ptotS(l,i+1,j,k)-ptotS(l,i,j,k))/dx
                    ro=half*(rS(l,i,j,k)+rS(l,i+1,j,k))
                    vo=half*(vS(l,i,j,k)+vS(l,i+1,j,k))
                    wo=half*(wS(l,i,j,k)+wS(l,i+1,j,k))
                    bo=half*(BS(l,i,j,k)+BS(l,i+1,j,k))
                    co=half*(CS(l,i,j,k)+CS(l,i+1,j,k))
                    pressure=(ro*(vo*vo+wo*wo)-bo*bo-co*co)*(ds(i+1,j,k,1)-ds(i,j,k,1))
                    flx_pre_tot(l,i,j,k)=flx_pre_tot(l,i,j,k) + pressure
                 end do
              end do
           end do
        end do
     endif
  endif

  if (xdim==2) then
     if (cartesian.or.cylindrical) then
        do k = klo, khi
           do j = jlo, jhi-1
              do i = ilo, ihi
                 do l = 1, ngrid
                    flx_pre_tot(l,i,j,k)=-(ptotS(l,i,j+1,k)*ds(i,j+1,k,2)-ptotS(l,i,j,k)*ds(i,j,k,2))
                 end do
              end do
           end do
        end do
     endif
     if (spherical) then
        do k = klo, khi
           do j = jlo, jhi-1
              cotanxc=cos(y(j))/sin(y(j))
              do i = ilo, ihi
                 do l = 1, ngrid
                    flx_pre_tot(l,i,j,k)=-dv(i,j,k)*(ptotS(l,i,j+1,k)-ptotS(l,i,j,k))/dy/x(i)
                    ro=half*(rS(l,i,j,k)+rS(l,i,j+1,k))
                    wo=half*(wS(l,i,j,k)+wS(l,i,j+1,k))
                    co=half*(cS(l,i,j,k)+cS(l,i,j+1,k))
                    pressure=(ro*wo*wo-co*co)*cotanxc*half*(ds(i+1,j,k,1)-ds(i,j,k,1))
                    flx_pre_tot(l,i,j,k)=flx_pre_tot(l,i,j,k) + pressure
                 end do
              end do
           end do
        end do
     endif
  endif

  if (xdim==3) then
     do k = klo, khi-1
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 flx_pre_tot(l,i,j,k)=-(ptotS(l,i,j,k+1)*ds(i,j,k+1,3)-ptotS(l,i,j,k)*ds(i,j,k,3))
              end do
           end do
        end do
     end do
  endif

  deallocate(qleft,qright,qgdnv,fgdnv,umean)
  deallocate(rS,uS,vS,wS,AS,BS,CS,ptotS)

  return

end subroutine cmpflxm

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CMP_MAG_FLX
!
!> 2D Riemann solver to compute EMF at cell edges.
!! Indices of the 2 planar velocity lp1 and lp2 and the
!! orthogonal one, lor and idem for the magnetic field.
!<
subroutine cmp_mag_flx(qRT,irt1,irt2,jrt1,jrt2,krt1,krt2, &
     &               qRB,irb1,irb2,jrb1,jrb2,krb1,krb2, &
     &               qLT,ilt1,ilt2,jlt1,jlt2,klt1,klt2, &
     &               qLB,ilb1,ilb2,jlb1,jlb2,klb1,klb2, &
     &                   ilo ,ihi ,jlo ,jhi ,klo ,khi , &
     &                   lp1 ,lp2 ,lor ,bp1 ,bp2 ,bor , &
     &               emf,ngrid)

  use amr_parameters
  use hydro_parameters
  use const
  use variables , only : x,cartesian,fargo

  implicit none

  integer :: ngrid
  integer :: lp1,lp2,lor,bp1,bp2,bor,nbvar
  integer :: irt1,irt2,jrt1,jrt2,krt1,krt2
  integer :: irb1,irb2,jrb1,jrb2,krb1,krb2
  integer :: ilt1,ilt2,jlt1,jlt2,klt1,klt2
  integer :: ilb1,ilb2,jlb1,jlb2,klb1,klb2
  integer :: ilo,ihi,jlo,jhi,klo,khi

  real(dp), dimension(1:nvector,irt1:irt2,jrt1:jrt2,krt1:krt2,1:nvar,1:3) :: qRT
  real(dp), dimension(1:nvector,irb1:irb2,jrb1:jrb2,krb1:krb2,1:nvar,1:3) :: qRB
  real(dp), dimension(1:nvector,ilt1:ilt2,jlt1:jlt2,klt1:klt2,1:nvar,1:3) :: qLT
  real(dp), dimension(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2,1:nvar,1:3) :: qLB
  real(dp), dimension(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2           ) :: emf

  ! local variables
  integer ::i,j,k,l,n,idim,xdim,m
  real(dp), dimension(:), allocatable :: qLL,qRL,qLR,qRR
  real(dp), dimension(:), allocatable :: qleft,qright,fmean_x,fmean_y,qtmp
  real(dp) :: ELL,ERL,ELR,ERR,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  real(dp) :: zero_flux,E,shear
  real(dp) :: cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  real(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  real(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  real(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  real(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR
  real(dp) :: ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  real(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  real(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  real(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  real(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  real(dp) :: AstarT,AstarB,BstarR,BstarL

  allocate(qLL(nvar),qRL(nvar),qLR(nvar),qRR(nvar))
  allocate(qleft(nvar),qright(nvar),fmean_x(nvar),fmean_y(nvar),qtmp(nvar))

  xdim = lor - 1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Density
              qLL(1) = qRT(l,i,j,k,1,xdim)
              qRL(1) = qLT(l,i,j,k,1,xdim)
              qLR(1) = qRB(l,i,j,k,1,xdim)
              qRR(1) = qLB(l,i,j,k,1,xdim)

              ! Pressure 
#if ISO==1
              qLL(2) = qLL(1)*ciso**2
              qRL(2) = qRL(1)*ciso**2
              qLR(2) = qLR(1)*ciso**2
              qRR(2) = qRR(1)*ciso**2
#else
              qLL(2) = qRT(l,i,j,k,5,xdim)
              qRL(2) = qLT(l,i,j,k,5,xdim)
              qLR(2) = qRB(l,i,j,k,5,xdim)
              qRR(2) = qLB(l,i,j,k,5,xdim)
#endif

              ! First parallel velocity 
              qLL(3) = qRT(l,i,j,k,lp1,xdim)
              qRL(3) = qLT(l,i,j,k,lp1,xdim)
              qLR(3) = qRB(l,i,j,k,lp1,xdim)
              qRR(3) = qLB(l,i,j,k,lp1,xdim)

              ! Second parallel velocity 
              qLL(4) = qRT(l,i,j,k,lp2,xdim)
              qRL(4) = qLT(l,i,j,k,lp2,xdim)
              qLR(4) = qRB(l,i,j,k,lp2,xdim)
              qRR(4) = qLB(l,i,j,k,lp2,xdim)

              ! First parallel magnetic field (enforce continuity)
              qLL(6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qRL(6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qLR(6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))
              qRR(6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))

              ! Second parallel magnetic field (enforce continuity)
              qLL(7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRL(7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))
              qLR(7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRR(7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))

              ! Orthogonal velocity 
              qLL(5) = qRT(l,i,j,k,lor,xdim)
              qRL(5) = qLT(l,i,j,k,lor,xdim)
              qLR(5) = qRB(l,i,j,k,lor,xdim)
              qRR(5) = qLB(l,i,j,k,lor,xdim)

              ! Orthogonal magnetic Field
              qLL(8) = qRT(l,i,j,k,bor,xdim)
              qRL(8) = qLT(l,i,j,k,bor,xdim)
              qLR(8) = qRB(l,i,j,k,bor,xdim)
              qRR(8) = qLB(l,i,j,k,bor,xdim)

              ! Compute final fluxes

              ! vx*by - vy*bx at the four edge centers
              ELL = qLL(3)*qLL(7) - qLL(4)*qLL(6)
              ERL = qRL(3)*qRL(7) - qRL(4)*qRL(6)
              ELR = qLR(3)*qLR(7) - qLR(4)*qLR(6)
              ERR = qRR(3)*qRR(7) - qRR(4)*qRR(6)

              if(iriemann2d==ihlld)then

                 rLL=qLL(1); pLL=qLL(2); uLL=qLL(3); vLL=qLL(4); ALL=qLL(6); BLL=qLL(7) ; CLL=qLL(8) 
                 rLR=qLR(1); pLR=qLR(2); uLR=qLR(3); vLR=qLR(4); ALR=qLR(6); BLR=qLR(7) ; CLR=qLR(8) 
                 rRL=qRL(1); pRL=qRL(2); uRL=qRL(3); vRL=qRL(4); ARL=qRL(6); BRL=qRL(7) ; CRL=qRL(8) 
                 rRR=qRR(1); pRR=qRR(2); uRR=qRR(3); vRR=qRR(4); ARR=qRR(6); BRR=qRR(7) ; CRR=qRR(8) 

                 ! Compute 4 fast magnetosonic velocity relative to x direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(3); qtmp(4)=qLL(6); qtmp(5)=qLL(4); qtmp(6)=qLL(7)
                 call find_speed_fast(qtmp,cfastLLx)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(3); qtmp(4)=qLR(6); qtmp(5)=qLR(4); qtmp(6)=qLR(7)
                 call find_speed_fast(qtmp,cfastLRx)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(3); qtmp(4)=qRL(6); qtmp(5)=qRL(4); qtmp(6)=qRL(7)
                 call find_speed_fast(qtmp,cfastRLx)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(3); qtmp(4)=qRR(6); qtmp(5)=qRR(4); qtmp(6)=qRR(7)
                 call find_speed_fast(qtmp,cfastRRx)

                 ! Compute 4 fast magnetosonic velocity relative to y direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(4); qtmp(4)=qLL(7); qtmp(5)=qLL(3); qtmp(6)=qLL(6)
                 call find_speed_fast(qtmp,cfastLLy)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(4); qtmp(4)=qLR(7); qtmp(5)=qLR(3); qtmp(6)=qLR(6)
                 call find_speed_fast(qtmp,cfastLRy)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(4); qtmp(4)=qRL(7); qtmp(5)=qRL(3); qtmp(6)=qRL(6)
                 call find_speed_fast(qtmp,cfastRLy)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(4); qtmp(4)=qRR(7); qtmp(5)=qRR(3); qtmp(6)=qRR(6)
                 call find_speed_fast(qtmp,cfastRRy)

                 SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                 SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                 SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
                 ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)

                 ELL=uLL*BLL-vLL*ALL
                 ELR=uLR*BLR-vLR*ALR
                 ERL=uRL*BRL-vRL*ARL
                 ERR=uRR*BRR-vRR*ARR

                 PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
                 PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
                 PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
                 PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

                 rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
                 rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
                 rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
                 rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

                 ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
                 vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

                 rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
                 rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
                 rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
                 EstarLLx=ustar*BstarLL-vLL  *ALL
                 EstarLLy=uLL  *BLL    -vstar*AstarLL
                 EstarLL =ustar*BstarLL-vstar*AstarLL

                 rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
                 rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
                 rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
                 EstarLRx=ustar*BstarLR-vLR  *ALR
                 EstarLRy=uLR  *BLR    -vstar*AstarLR
                 EstarLR =ustar*BstarLR-vstar*AstarLR

                 rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
                 rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
                 rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
                 EstarRLx=ustar*BstarRL-vRL  *ARL
                 EstarRLy=uRL  *BRL    -vstar*AstarRL
                 EstarRL =ustar*BstarRL-vstar*AstarRL

                 rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
                 rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
                 rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
                 EstarRRx=ustar*BstarRR-vRR  *ARR
                 EstarRRy=uRR  *BRR    -vstar*AstarRR
                 EstarRR =ustar*BstarRR-vstar*AstarRR

                 calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                      &       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
                 calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                      &       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
                 calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                      &       abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
                 calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                      &       abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
                 SAL=min(ustar-calfvenL,zero); SAR=max(ustar+calfvenR,zero)
                 SAB=min(vstar-calfvenB,zero); SAT=max(vstar+calfvenT,zero)
                 AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
                 BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

                 if(SB>0d0)then
                    if(SL>0d0)then
                       E=ELL
                    else if(SR<0d0)then
                       E=ERL
                    else
                       E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
                    endif
                 else if (ST<0d0)then
                    if(SL>0d0)then
                       E=ELR
                    else if(SR<0d0)then
                       E=ERR
                    else
                       E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
                    endif
                 else if(SL>0d0)then
                    E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
                 else if(SR<0d0)then
                    E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
                 else
                    E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
                         & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
                 endif

                 emf(l,i,j,k) = E

              else if(iriemann2d==ihllf)then

                 ! Compute 4 fast magnetosonic velocity relative to x direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(3); qtmp(4)=qLL(6); qtmp(5)=qLL(4); qtmp(6)=qLL(7)
                 vLLx=qtmp(3); call find_speed_fast(qtmp,cLLx)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(3); qtmp(4)=qLR(6); qtmp(5)=qLR(4); qtmp(6)=qLR(7)
                 vLRx=qtmp(3); call find_speed_fast(qtmp,cLRx)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(3); qtmp(4)=qRL(6); qtmp(5)=qRL(4); qtmp(6)=qRL(7)
                 vRLx=qtmp(3); call find_speed_fast(qtmp,cRLx)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(3); qtmp(4)=qRR(6); qtmp(5)=qRR(4); qtmp(6)=qRR(7)
                 vRRx=qtmp(3); call find_speed_fast(qtmp,cRRx)

                 ! Compute 4 fast magnetosonic velocity relative to y direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(4); qtmp(4)=qLL(7); qtmp(5)=qLL(3); qtmp(6)=qLL(6)
                 vLLy=qtmp(3); call find_speed_fast(qtmp,cLLy)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(4); qtmp(4)=qLR(7); qtmp(5)=qLR(3); qtmp(6)=qLR(6)
                 vLRy=qtmp(3); call find_speed_fast(qtmp,cLRy)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(4); qtmp(4)=qRL(7); qtmp(5)=qRL(3); qtmp(6)=qRL(6)
                 vRLy=qtmp(3); call find_speed_fast(qtmp,cRLy)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(4); qtmp(4)=qRR(7); qtmp(5)=qRR(3); qtmp(6)=qRR(6)
                 vRRy=qtmp(3); call find_speed_fast(qtmp,cRRy)

                 SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx,smallc),zero)
                 SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx,smallc),zero)
                 SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy,smallc),zero)
                 ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy,smallc),zero)

                 emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                      -ST*SB/(ST-SB)*(qRR(6)-qLL(6)) &
                      +SR*SL/(SR-SL)*(qRR(7)-qLL(7))

              else if (iriemann2d==ihlla)then

                 ! Compute 4 Alfven velocity relative to x direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(3); qtmp(4)=qLL(6); qtmp(5)=qLL(4); qtmp(6)=qLL(7)
                 vLLx=qtmp(3); call find_speed_alfven(qtmp,cLLx)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(3); qtmp(4)=qLR(6); qtmp(5)=qLR(4); qtmp(6)=qLR(7)
                 vLRx=qtmp(3); call find_speed_alfven(qtmp,cLRx)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(3); qtmp(4)=qRL(6); qtmp(5)=qRL(4); qtmp(6)=qRL(7)
                 vRLx=qtmp(3); call find_speed_alfven(qtmp,cRLx)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(3); qtmp(4)=qRR(6); qtmp(5)=qRR(4); qtmp(6)=qRR(7)
                 vRRx=qtmp(3); call find_speed_alfven(qtmp,cRRx)
                 ! Compute 4 Alfven relative to y direction
                 qtmp(1)=qLL(1); qtmp(2)=qLL(2); qtmp(7)=qLL(5); qtmp(8)=qLL(8)
                 qtmp(3)=qLL(4); qtmp(4)=qLL(7); qtmp(5)=qLL(3); qtmp(6)=qLL(6)
                 vLLy=qtmp(3); call find_speed_alfven(qtmp,cLLy)
                 qtmp(1)=qLR(1); qtmp(2)=qLR(2); qtmp(7)=qLR(5); qtmp(8)=qLR(8)
                 qtmp(3)=qLR(4); qtmp(4)=qLR(7); qtmp(5)=qLR(3); qtmp(6)=qLR(6)
                 vLRy=qtmp(3); call find_speed_alfven(qtmp,cLRy)
                 qtmp(1)=qRL(1); qtmp(2)=qRL(2); qtmp(7)=qRL(5); qtmp(8)=qRL(8)
                 qtmp(3)=qRL(4); qtmp(4)=qRL(7); qtmp(5)=qRL(3); qtmp(6)=qRL(6)
                 vRLy=qtmp(3); call find_speed_alfven(qtmp,cRLy)
                 qtmp(1)=qRR(1); qtmp(2)=qRR(2); qtmp(7)=qRR(5); qtmp(8)=qRR(8)
                 qtmp(3)=qRR(4); qtmp(4)=qRR(7); qtmp(5)=qRR(3); qtmp(6)=qRR(6)
                 vRRy=qtmp(3); call find_speed_alfven(qtmp,cRRy)

                 SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx,smallc),zero)
                 SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx,smallc),zero)
                 SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy,smallc),zero)
                 ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy,smallc),zero)

                 emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                      -ST*SB/(ST-SB)*(qRR(6)-qLL(6)) &
                      +SR*SL/(SR-SL)*(qRR(7)-qLL(7))

              else

                 ! find the average value of E
                 E = forth*(ELL+ERL+ELR+ERR)

                 ! call the first solver in the x direction
                 ! density
                 qleft (1) = half*(qLL(1)+qLR(1))
                 qright(1) = half*(qRR(1)+qRL(1))

                 ! pressure
                 qleft (2) = half*(qLL(2)+qLR(2))
                 qright(2) = half*(qRR(2)+qRL(2))

                 ! vt1 becomes normal velocity
                 qleft (3) = half*(qLL(3)+qLR(3))
                 qright(3) = half*(qRR(3)+qRL(3))

                 ! bt1 becomes normal magnetic field
                 qleft (4) = half*(qLL(6)+qLR(6))
                 qright(4) = half*(qRR(6)+qRL(6))

                 ! vt2 becomes transverse velocity field
                 qleft (5) = half*(qLL(4)+qLR(4))
                 qright(5) = half*(qRR(4)+qRL(4))

                 ! bt2 becomes transverse magnetic field 
                 qleft (6) = half*(qLL(7)+qLR(7))
                 qright(6) = half*(qRR(7)+qRL(7))

                 ! velocity component perp. to the plane is now transverse
                 qleft (7) = half*(qLL(5)+qLR(5))
                 qright(7) = half*(qRR(5)+qRL(5))

                 ! magnetic field component perp. to the plane is now transverse
                 qleft (8) = half*(qLL(8)+qLR(8))
                 qright(8) = half*(qRR(8)+qRL(8))

                 zero_flux = 0.0
                 select case (iriemann2d)
                 case (iroe)
                    call athena_roe   (qleft,qright,fmean_x,zero_flux)
                 case (illf)
                    call lax_friedrich(qleft,qright,fmean_x,zero_flux)
                 case (iupwind)
                    call upwind       (qleft,qright,fmean_x,zero_flux)
                 end select

                 ! call the second solver in the y direction
                 ! density
                 qleft (1) = half*(qLL(1)+qRL(1))
                 qright(1) = half*(qRR(1)+qLR(1))

                 ! pressure
                 qleft (2) = half*(qLL(2)+qRL(2))
                 qright(2) = half*(qRR(2)+qLR(2))

                 ! vt2 becomes normal velocity
                 qleft (3) = half*(qLL(4)+qRL(4))
                 qright(3) = half*(qRR(4)+qLR(4))

                 ! bt2 becomes normal magnetic field
                 qleft (4) = half*(qLL(7)+qRL(7))
                 qright(4) = half*(qRR(7)+qLR(7))

                 ! vt1 becomes transverse velocity field 
                 qleft (5) = half*(qLL(3)+qRL(3))
                 qright(5) = half*(qRR(3)+qLR(3))

                 ! bt1 becomes transverse magnetic field 
                 qleft (6) = half*(qLL(6)+qRL(6))
                 qright(6) = half*(qRR(6)+qLR(6))

                 ! velocity component perp. to the plane is now transverse
                 qleft (7) = half*(qLL(5)+qRL(5))
                 qright(7) = half*(qRR(5)+qLR(5))

                 ! magnetic field component perp. to the plane is now transverse
                 qleft (8) = half*(qLL(8)+qRL(8))
                 qright(8) = half*(qRR(8)+qLR(8))

                 zero_flux = 0.
                 select case (iriemann2d)
                 case (iroe)
                    call athena_roe   (qleft,qright,fmean_y,zero_flux)
                 case (illf)
                    call lax_friedrich(qleft,qright,fmean_y,zero_flux)
                 case (iupwind)
                    call upwind       (qleft,qright,fmean_y,zero_flux)
                 end select

                 ! compute the final value of E including the 2D diffusive
                 ! terms that ensure stability 
                 emf(l,i,j,k) = E + (fmean_x(6) - fmean_y(6))

              endif

              ! Upwind solver in case of the shearing box
              if ((Omega0>0.d0).and.(cartesian).and..not.(fargo)) then
                 if (xdim==1) then
                    shear=-1.5*Omega0*x(i)
                    if (shear>0.d0) then
                       emf(l,i,j,k)=emf(l,i,j,k) + shear*qLL(7)
                    else
                       emf(l,i,j,k)=emf(l,i,j,k) + shear*qRR(7)
                    endif
                 endif
                 if (xdim==3) then
                    shear=-1.5*Omega0*half*(x(i)+x(i-1))
                    if (shear>0.d0) then
                       emf(l,i,j,k)=emf(l,i,j,k) - shear*qLL(6)
                    else
                       emf(l,i,j,k)=emf(l,i,j,k) - shear*qRR(6)
                    endif
                 endif
              endif

           end do
        end do
     end do
  end do

  deallocate(qLL,qRL,qLR,qRR)
  deallocate(qleft,qright,fmean_x,fmean_y,qtmp)

  return

end subroutine cmp_mag_flx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!  Subroutine USLOPE
!
!> Computes mhd slopes.
!<
subroutine uslope(bf,q,dq,dbf,ngrid)

  use amr_parameters
  use hydro_parameters
  use const

  implicit none

  integer :: ngrid
  real(dp), dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3          ) :: bf
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar       ) :: q 
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: dq
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:3   ,1:ndim) :: dbf

  ! local variables
  integer :: i,j,k,l,n
  integer :: ilo,ihi,jlo,jhi,klo,khi,ind_n
  real(dp)    :: dsgn,dlim,dcen,dlft,drgt,slop
  real(dp)    :: dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)    :: dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)    :: dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)    :: dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)    :: vmin,vmax,dfx,dfy,dfz,dff,xslope_type

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     dbf=zero
     return
  end if

  !1D
  if (ndim==1) then
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    if(slope_type==1.or.slope_type==2)then  ! minmod or average
                       dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                       drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one, dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                    else
                       write(*,*)'Unknown slope type'
                       stop
                    end if
                 end do
              end do
           end do
        end do
     end do
  endif

  !2D
  if (ndim==2) then
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        do n = 1, nvar
           do k = klo, khi
              do j = jlo, jhi
                 do i = ilo, ihi
                    do l = 1, ngrid
                       ! slopes in first coordinate direction
                       dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                       drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one, dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                       ! slopes in second coordinate direction
                       dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                       drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one,dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                    end do
                 end do
              end do
           end do
        end do

     else if(slope_type==3)then ! positivity preserving 2d unsplit slope

        do n = 1, nvar
           do k = klo, khi
              do j = jlo, jhi
                 do i = ilo, ihi
                    do l= 1, ngrid

                       dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                       dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                       dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                       dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                       dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                       dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                       dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                       dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                       dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)

                       vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                       vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

                       dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                       dff  = half*(abs(dfx)+abs(dfy))

                       if(dff>zero)then
                          slop = min(one,min(abs(vmin),abs(vmax))/dff)
                       else
                          slop = one
                       endif

                       dlim = slop

                       dq(l,i,j,k,n,1) = dlim*dfx
                       dq(l,i,j,k,n,2) = dlim*dfy

                    end do
                 end do
              end do
           end do
        end do
     else
        write(*,*)"Unknown slope type"
        stop
     endif

     ! 1D transverse TVD slopes for face-centered magnetic fields
     ! Bx along direction Y 
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1 ! WARNING HERE
              do l = 1, ngrid
                 dlft = slope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
                 drgt = slope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,1,2) = dsgn*min(dlim,abs(dcen))
              end do
           enddo
        end do
     end do
     ! By along direction X
     do k = klo, khi
        do j = jlo, jhi+1 ! WARNING HERE
           do i = ilo, ihi
              do l = 1, ngrid
                 dlft = slope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
                 drgt = slope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
              end do
           enddo
        end do
     end do
  endif

  !3D
  if (ndim==3) then
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        do n = 1, nvar
           do k = klo, khi
              do j = jlo, jhi
                 do i = ilo, ihi
                    do l = 1, ngrid
                       ! slopes in first coordinate direction
                       dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                       drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one, dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                       ! slopes in second coordinate direction
                       dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                       drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one,dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                       ! slopes in third coordinate direction
                       dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                       drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                       dcen = half*(dlft+drgt)/slope_type
                       dsgn = sign(one,dcen)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                    end do
                 end do
              end do
           end do
        end do

     else if(slope_type==3)then ! positivity preserving 3d unsplit slope
        do n = 1, nvar
           do k = klo, khi
              do j = jlo, jhi
                 do i = ilo, ihi
                    do l = 1, ngrid
                       dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                       dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                       dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                       dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                       dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                       dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                       dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                       dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                       dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                       dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                       dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                       dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                       dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                       dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                       dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                       dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                       dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                       dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                       dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                       dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                       dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                       dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                       dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                       dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                       dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                       dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                       dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                       vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                            &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                            &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                       vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                            &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                            &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                       dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                       dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                       dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

                       if(dff>zero)then
                          slop = min(one,min(abs(vmin),abs(vmax))/dff)
                       else
                          slop = one
                       endif

                       dlim = slop

                       dq(l,i,j,k,n,1) = dlim*dfx
                       dq(l,i,j,k,n,2) = dlim*dfy
                       dq(l,i,j,k,n,3) = dlim*dfz

                    end do
                 end do
              end do
           end do
        end do
     else
        write(*,*)'Unknown slope type'
        stop
     endif

     ! 2D transverse TVD slopes for face-centered magnetic fields
     ! Bx along direction Y and Z
     xslope_type=min(slope_type,2)
     if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi+1 ! WARNING HERE
                 do l = 1, ngrid
                    ! slopes along Y
                    dlft = xslope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
                    drgt = xslope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,1,2) = dsgn*min(dlim,abs(dcen))
                    ! slopes along Z
                    dlft = xslope_type*(bf(l,i,j,k  ,1) - bf(l,i,j,k-1,1))
                    drgt = xslope_type*(bf(l,i,j,k+1,1) - bf(l,i,j,k  ,1))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,1,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     endif
     ! By along direction X and Z
     if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
        do k = klo, khi
           do j = jlo, jhi+1 ! WARNING HERE
              do i = ilo, ihi
                 do l = 1, ngrid
                    ! slopes along X
                    dlft = xslope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
                    drgt = xslope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
                    ! slopes along Z
                    dlft = xslope_type*(bf(l,i,j,k  ,2) - bf(l,i,j,k-1,2))
                    drgt = xslope_type*(bf(l,i,j,k+1,2) - bf(l,i,j,k  ,2))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,2,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     endif
     ! Bz along direction X and Y
     if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
        do k = klo, khi+1 ! WARNING HERE
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    ! slopes along X
                    dlft = xslope_type*(bf(l,i  ,j,k,3) - bf(l,i-1,j,k,3))
                    drgt = xslope_type*(bf(l,i+1,j,k,3) - bf(l,i  ,j,k,3))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,3,1) = dsgn*min(dlim,abs(dcen))
                    ! slopes along Y
                    dlft = xslope_type*(bf(l,i,j  ,k,3) - bf(l,i,j-1,k,3))
                    drgt = xslope_type*(bf(l,i,j+1,k,3) - bf(l,i,j  ,k,3))
                    dcen = half*(dlft+drgt)/xslope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dbf(l,i,j,k,3,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     endif
  endif

  return
end subroutine uslope
