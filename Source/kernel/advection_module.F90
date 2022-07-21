module advection_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: compute_flux

contains

  subroutine compute_flux(q, qd_lo, qd_hi, lo, hi, dx, flux1, flux2, flux3)

    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    use nc_module, only : urho, umx, umy, umz, ueden, nvar, &
      qrho,qu,qv,qw,qp,qvar, smallp, smallr, gamma
    use fluxsplit_module, only : flux_split

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(  out) :: flux1(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,qvar)
    real(rt), intent(  out) :: flux2(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,qvar)
    real(rt), intent(  out) :: flux3(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,qvar)

    real(rt), dimension(5) :: qL, qR, flux_loc
    real(rt) :: tmp
    integer :: i,j,k

    ! X-direction
    do k = lo(3) ,hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          qL = q(i-1,j,k, :) + &
            0.5d0 * minmod(q(i,j,k,:)-q(i-1,j,k,:), q(i-1,j,k,:)-q(i-2,j,k,:))
          qR = q(i,j,k, :) - &
            0.5d0 * minmod(q(i,j,k,:)-q(i-1,j,k,:), q(i+1,j,k,:)-q(i,j,k,:))

          call flux_split(qL, qR, flux_loc)

          flux1(i,j,k,:) = flux_loc

        end do
      end do
    end do

    ! Y-direction
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
          qL = q(i,j-1,k, :) + &
            0.5d0 * minmod(q(i,j,k,:)-q(i,j-1,k,:), q(i,j-1,k,:)-q(i,j-2,k,:))
          qR = q(i,j,k, :) - &
            0.5d0 * minmod(q(i,j,k,:)-q(i,j-1,k,:), q(i,j+1,k,:)-q(i,j,k,:))

          ! rotate here
          tmp = qL(qu)
          qL(qu) = qL(qv)
          qL(qv) = -tmp

          tmp = qR(qu)
          qR(qu) = qR(qv)
          qR(qv) = -tmp
          call flux_split(qL, qR, flux_loc)

          ! revert
          tmp = flux_loc(qu)
          flux_loc(qu) = -flux_loc(qv)
          flux_loc(qv) = tmp

          flux2(i,j,k,:) = flux_loc

        end do
      end do
    end do

    ! Z-direction
    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          qL = q(i,j,k-1, :) + &
            0.5d0 * minmod(q(i,j,k,:)-q(i,j,k-1,:), q(i,j,k-1,:)-q(i,j,k-2,:))
          qR = q(i,j,k, :) - &
            0.5d0 * minmod(q(i,j,k,:)-q(i,j,k-1,:), q(i,j,k+1,:)-q(i,j,k,:))
            
          ! exchange u and w
          tmp = qL(qu)
          qL(qu) = qL(qw)
          qL(qw) = -tmp

          tmp = qR(qu)
          qR(qu) = qR(qw)
          qR(qw) = -tmp

          call flux_split(qL, qR, flux_loc)

          ! revert
          tmp = flux_loc(qu)
          flux_loc(qu) = -flux_loc(qw)
          flux_loc(qw) = tmp

          flux3(i,j,k,:) = flux_loc

        end do
      end do
    end do

  end subroutine compute_flux

  pure function minmod(a, b) result(res)
    real(rt) , intent(in), dimension(5) :: a, b
    real(rt), dimension(5) :: res
    if (a(1)*b(1) > 0) then
      if (abs(a(1))>abs(b(1))) then
        res = b
      else
        res = a
      end if
    else
      res = 0.d0
    endif
  end function minmod

end module advection_module