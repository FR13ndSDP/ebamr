module viscous_module
  use amrex_fort_module, only : rt=>amrex_real
  use nc_module, only : gamma, cv, R, pr
  implicit none
  private
  public :: compute_viscous_flux

contains

  subroutine compute_viscous_flux(q, qd_lo, qd_hi, lo, hi, dx, flux1, flux2, flux3)
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    use nc_module, only : urho, umx, umy, umz, ueden, nvar, &
      qrho,qu,qv,qw,qp,qvar, smallp, smallr
    
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(  out) :: flux1(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,qvar)
    real(rt), intent(  out) :: flux2(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,qvar)
    real(rt), intent(  out) :: flux3(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,qvar)

    integer :: i,j,k
    real(rt) :: dxinv(3), Prinv = 1.d0/Pr
    real(rt) :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, muf
    real(rt) :: dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz, divu, dTdx
    real(rt) :: dvdz, dwdy, dTdy, dTdz
    real(rt), parameter :: twoThirds = 2.d0/3.d0
    real(rt) :: mu_0 = 1.716d-5, T_ref = 273.15d0, T_0 = 111.d0

    real(rt), dimension(:,:,:), pointer, contiguous :: T, mut, kt
    call amrex_allocate(T,   qd_lo(1), qd_hi(1), qd_lo(2), qd_hi(2), qd_lo(3),qd_hi(3))
    call amrex_allocate(mut, qd_lo(1), qd_hi(1), qd_lo(2), qd_hi(2), qd_lo(3),qd_hi(3))
    call amrex_allocate(kt,  qd_lo(1), qd_hi(1), qd_lo(2), qd_hi(2), qd_lo(3),qd_hi(3))

    dxinv = 1.d0/dx

    ! first compute temperature
    T = q(:,:,:, qp)/(R * q(:,:,:,qrho))
    ! derive mu
    mut = mu_0*(T/T_ref)**1.5d0*(T_ref+T_0)/(T+T_0)
    ! derive k
    kt = mut * cv * gamma * Prinv

    ! x-direction
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             dTdx = (T(i,j,k)-T(i-1,j,k))*dxinv(1)
             dudx = (q(i,j,k,qu)-q(i-1,j,k,qu))*dxinv(1)
             dvdx = (q(i,j,k,qv)-q(i-1,j,k,qv))*dxinv(1)
             dwdx = (q(i,j,k,qw)-q(i-1,j,k,qw))*dxinv(1)
             dudy = (q(i,j+1,k,qu)+q(i-1,j+1,k,qu)-q(i,j-1,k,qu)-q(i-1,j-1,k,qu))*(0.25d0*dxinv(2))
             dvdy = (q(i,j+1,k,qv)+q(i-1,j+1,k,qv)-q(i,j-1,k,qv)-q(i-1,j-1,k,qv))*(0.25d0*dxinv(2))
             dudz = (q(i,j,k+1,qu)+q(i-1,j,k+1,qu)-q(i,j,k-1,qu)-q(i-1,j,k-1,qu))*(0.25d0*dxinv(3))
             dwdz = (q(i,j,k+1,qw)+q(i-1,j,k+1,qw)-q(i,j,k-1,qw)-q(i-1,j,k-1,qw))*(0.25d0*dxinv(3))
             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mut(i,j,k)+mut(i-1,j,k))
             tauxx = muf*(2.d0*dudx-twoThirds*divu)
             tauxy = muf*(dudy+dvdx)
             tauxz = muf*(dudz+dwdx)
             flux1(i,j,k,urho)  = 0.d0
             flux1(i,j,k,umx)   = -tauxx
             flux1(i,j,k,umy)   = -tauxy
             flux1(i,j,k,umz)   = -tauxz
             flux1(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i-1,j,k,qu))*tauxx &
                  &                      +(q(i,j,k,qv)+q(i-1,j,k,qv))*tauxy &
                  &                      +(q(i,j,k,qw)+q(i-1,j,k,qw))*tauxz &
                  &                      +(kt(i,j,k) +kt(i-1,j,k))*dTdx)
          end do
       end do
    end do

    ! y-direction
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             dTdy = (T(i,j,k)-T(i,j-1,k))*dxinv(2)
             dudy = (q(i,j,k,qu)-q(i,j-1,k,qu))*dxinv(2)
             dvdy = (q(i,j,k,qv)-q(i,j-1,k,qv))*dxinv(2)
             dwdy = (q(i,j,k,qw)-q(i,j-1,k,qw))*dxinv(2)
             dudx = (q(i+1,j,k,qu)+q(i+1,j-1,k,qu)-q(i-1,j,k,qu)-q(i-1,j-1,k,qu))*(0.25d0*dxinv(1))
             dvdx = (q(i+1,j,k,qv)+q(i+1,j-1,k,qv)-q(i-1,j,k,qv)-q(i-1,j-1,k,qv))*(0.25d0*dxinv(1))
             dvdz = (q(i,j,k+1,qv)+q(i,j-1,k+1,qv)-q(i,j,k-1,qv)-q(i,j-1,k-1,qv))*(0.25d0*dxinv(3))
             dwdz = (q(i,j,k+1,qw)+q(i,j-1,k+1,qw)-q(i,j,k-1,qw)-q(i,j-1,k-1,qw))*(0.25d0*dxinv(3))
             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mut(i,j,k)+mut(i,j-1,k))
             tauyy = muf*(2.d0*dvdy-twoThirds*divu)
             tauxy = muf*(dudy+dvdx)
             tauyz = muf*(dwdy+dvdz)
             flux2(i,j,k,urho)  = 0.d0
             flux2(i,j,k,umx)   = -tauxy
             flux2(i,j,k,umy)   = -tauyy
             flux2(i,j,k,umz)   = -tauyz
             flux2(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i,j-1,k,qu))*tauxy &
                  &                      +(q(i,j,k,qv)+q(i,j-1,k,qv))*tauyy &
                  &                      +(q(i,j,k,qw)+q(i,j-1,k,qw))*tauyz &
                  &                      +(kt(i,j,k) +kt(i,j-1,k))*dTdy)
          end do
       end do
    end do

    ! z-direction
    do       k = lo(3), hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dTdz = (T(i,j,k)-T(i,j,k-1))*dxinv(3)
             dudz = (q(i,j,k,qu)-q(i,j,k-1,qu))*dxinv(3)
             dvdz = (q(i,j,k,qv)-q(i,j,k-1,qv))*dxinv(3)
             dwdz = (q(i,j,k,qw)-q(i,j,k-1,qw))*dxinv(3)
             dudx = (q(i+1,j,k,qu)+q(i+1,j,k-1,qu)-q(i-1,j,k,qu)-q(i-1,j,k-1,qu))*(0.25d0*dxinv(1))
             dwdx = (q(i+1,j,k,qw)+q(i+1,j,k-1,qw)-q(i-1,j,k,qw)-q(i-1,j,k-1,qw))*(0.25d0*dxinv(1))
             dvdy = (q(i,j+1,k,qv)+q(i,j+1,k-1,qv)-q(i,j-1,k,qv)-q(i,j-1,k-1,qv))*(0.25d0*dxinv(2))
             dwdy = (q(i,j+1,k,qw)+q(i,j+1,k-1,qw)-q(i,j-1,k,qw)-q(i,j-1,k-1,qw))*(0.25d0*dxinv(2))
             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mut(i,j,k)+mut(i,j,k-1))
             tauxz = muf*(dudz+dwdx)
             tauyz = muf*(dvdz+dwdy)
             tauzz = muf*(2.d0*dwdz-twoThirds*divu)
             flux3(i,j,k,urho)  = 0.d0
             flux3(i,j,k,umx)   = -tauxz
             flux3(i,j,k,umy)   = -tauyz
             flux3(i,j,k,umz)   = -tauzz
             flux3(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i,j,k-1,qu))*tauxz &
                  &                      +(q(i,j,k,qv)+q(i,j,k-1,qv))*tauyz &
                  &                      +(q(i,j,k,qw)+q(i,j,k-1,qw))*tauzz &
                  &                      +(kt(i,j,k) +kt(i,j,k-1))*dTdz)
          end do
       end do
    end do

    call amrex_deallocate(T)
    call amrex_deallocate(mut)
    call amrex_deallocate(kt)
  end subroutine compute_viscous_flux
end module viscous_module
