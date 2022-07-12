module bc_fill_module

  implicit none

  public

contains

  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.

  ! subroutine nc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
  !   bind(C, name="nc_hypfill")

  !   use nc_module, only: NVAR
  !   use amrex_fort_module, only: dim=>amrex_spacedim
  !   use amrex_filcc_module
  !   use amrex_bc_types_module, only : amrex_bc_ext_dir
  !   use probdata_module
  !   use nc_physics_module, only: gamma, cv

  !   implicit none

  !   integer          :: adv_lo(3),adv_hi(3)
  !   integer          :: bc(dim,2,*)
  !   integer          :: domlo(3), domhi(3)
  !   double precision :: delta(3), xlo(3), time
  !   double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

  !   integer          :: i,j,k,n

  !   do n = 1,NVAR
  !     call amrex_filcc(adv(:,:,:,n),adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3),domlo,domhi,delta,xlo,bc(:,:,n))
  !   enddo

  !   !bc (direction, lo_or_hi, var)
  !   if ( bc(1,1,1).eq.amrex_bc_ext_dir .and. adv_lo(1).lt.domlo(1)) then
  !     do       k = adv_lo(3), adv_hi(3)
  !       do    j = adv_lo(2), adv_hi(2)
  !         do i = adv_lo(1), domlo(1)
  !           adv(i,j,k,1) = rho_l
  !           adv(i,j,k,2) = rho_l * u_l
  !           adv(i,j,k,3) = 0.d0
  !           adv(i,j,k,4) = 0.d0
  !           adv(i,j,k,5) = p_l/(gamma-1.d0) + 0.5d0*u_l*u_l*rho_l
  !           adv(i,j,k,6) = p_l/(gamma-1.d0)
  !           adv(i,j,k,7) = adv(i,j,k,6)/(rho_l * cv)
  !         end do
  !       end do
  !     end do
  !   end if

  ! end subroutine nc_hypfill



  ! subroutine nc_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
  !   bind(C, name="nc_denfill")

  !   use amrex_fort_module, only: dim=>amrex_spacedim
  !   use amrex_filcc_module
  !   use amrex_bc_types_module, only : amrex_bc_ext_dir
  !   use probdata_module
  !   implicit none

  !   include 'AMReX_bc_types.fi'

  !   integer          :: adv_lo(3),adv_hi(3)
  !   integer          :: bc(dim,2)
  !   integer          :: domlo(3), domhi(3)
  !   double precision :: delta(3), xlo(3), time
  !   double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

  !   integer :: i,j,k
  !   call amrex_filcc(adv,adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3),domlo,domhi,delta,xlo,bc)

  !   if ( bc(1,1).eq.amrex_bc_ext_dir .and. adv_lo(1).lt.domlo(1)) then
  !     do       k = adv_lo(3), adv_hi(3)
  !       do    j = adv_lo(2), adv_hi(2)
  !         do i = adv_lo(1), domlo(1)
  !           adv(i,j,k) = rho_l
  !         end do
  !       end do
  !     end do
  !   end if

  ! end subroutine nc_denfill

  subroutine nullfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
    bind(C, name="nullfill")
    use amrex_fort_module, only: dim=>amrex_spacedim
    use amrex_error_module, only : amrex_error
    implicit none
    include 'AMReX_bc_types.fi'
    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    call amrex_error("How did this happen?")
  end subroutine nullfill

end module bc_fill_module
