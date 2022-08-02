module advection_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: compute_flux
  ! TODO: verify weno_5th
contains

  subroutine compute_flux(q, qd_lo, qd_hi, lo, hi, dx, flux1, flux2, flux3)
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    use nc_module, only : urho, umx, umy, umz, ueden, nvar, &
      qrho,qu,qv,qw,qp,qvar, smallp, smallr, gamma, refactor_scheme
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
    integer :: i,j,k,n

    ! X-direction
    do k = lo(3) ,hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          if (refactor_scheme .eq. 1) then
            do n = 1,QVAR
              qL(n) = q(i-1,j,k, n) + &
                0.5d0 * minmod(q(i,j,k,n)-q(i-1,j,k,n), q(i-1,j,k,n)-q(i-2,j,k,n))
              qR(n) = q(i,j,k,n) - &
                0.5d0 * minmod(q(i,j,k,n)-q(i-1,j,k,n), q(i+1,j,k,n)-q(i,j,k,n))
            end do
          else
            call weno_5th(qL, qR, q(i-3:i+2,j,k,:))
          end if
          call flux_split(qL, qR, flux_loc)

          flux1(i,j,k,:) = flux_loc

        end do
      end do
    end do

    ! Y-direction
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
          if (refactor_scheme .eq. 1) then
            do n = 1,QVAR
              qL(n) = q(i,j-1,k, n) + &
                0.5d0 * minmod(q(i,j,k,n)-q(i,j-1,k,n), q(i,j-1,k,n)-q(i,j-2,k,n))
              qR(n) = q(i,j,k, n) - &
                0.5d0 * minmod(q(i,j,k,n)-q(i,j-1,k,n), q(i,j+1,k,n)-q(i,j,k,n))
            end do
          else 
            call weno_5th(qL, qR, q(i,j-3:j+2,k,:))
          end if
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
          if (refactor_scheme .eq. 1) then
            do n = 1,QVAR
              qL(n) = q(i,j,k-1, n) + &
                0.5d0 * minmod(q(i,j,k,n)-q(i,j,k-1,n), q(i,j,k-1,n)-q(i,j,k-2,n))
              qR(n) = q(i,j,k, n) - &
                0.5d0 * minmod(q(i,j,k,n)-q(i,j,k-1,n), q(i,j,k+1,n)-q(i,j,k,n))
            end do
          else
            call weno_5th(qL, qR, q(i,j,k-3:k+2,:))
          end if
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

  subroutine weno_5th(qL, qR, q)
    use nc_module, only : QVAR
    real(rt), intent(inout) :: qL(qvar), qR(qvar)
    real(rt), intent(in) :: q(-3:2, qvar)
    real(rt) :: flux(3, qvar), IS(3, qvar), alpha(3, qvar), omega(3, qvar), c(3)
    real(rt), parameter :: eps = 1d-6, p = 2.d0

    real(rt), parameter :: onethree = 1.d0/3.d0, sevensixth = 7.d0/6.d0, &
                           elvsixth = 11.d0/6.d0, onesixth = 1.d0/6.d0, &
                           fivesixth = 5.d0/6.d0, onefour = 1.d0/4.d0, &
                           thirtwelve  = 13.d0/12.d0

    c = (/1.d0/10.d0, 3.d0/5.d0, 3.d0/10.d0/)
    ! qL
    flux(1, :) = onethree * q(-3, :) - sevensixth*q(-2, :)+elvsixth*q(-1, :)
    flux(2, :) = -onesixth*q(-2,:) + fivesixth*q(-1, :) + onethree*q(0, :)
    flux(3,:) = onethree*q(-1,:)+fivesixth*q(0,:)-onesixth*q(1, :)

    IS(1, :) = onefour*(q(-3,:)-4.d0*q(-2,:)+3.d0*q(-1,:))**2 + &
               thirtwelve*(q(-3,:)-2.d0*q(-2,:)+q(-1,:))**2

    IS(2, :) = onefour*(q(-2,:)-q(0,:))**2 + &
               thirtwelve*(q(-2,:)-2.d0*q(-1,:)+q(0,:))**2
               
    IS(3, :) = onefour*(3.d0*q(-1,:) - 4.d0*q(0,:) +q(1, :))**2 + &
               thirtwelve*(q(-1,:)-2.d0*q(0,:)+q(1,:))**2

    alpha(1, :) = c(1)/(eps + IS(1,:))**p
    alpha(2, :) = c(2)/(eps + IS(2,:))**p
    alpha(3, :) = c(3)/(eps + IS(3,:))**p

    omega(1,:) = alpha(1, :)/(sum(alpha, dim=1))
    omega(2,:) = alpha(2, :)/(sum(alpha, dim=1))
    omega(3,:) = alpha(3, :)/(sum(alpha, dim=1))

    qL = omega(1,:)*flux(1,:) + omega(2,:)*flux(2,:) + omega(3,:)*flux(3,:)

    ! qR
    flux(1, :) = onethree * q(2, :) - sevensixth*q(1, :)+elvsixth*q(0, :)
    flux(2, :) = -onesixth*q(1,:) + fivesixth*q(0, :) + onethree*q(-1, :)
    flux(3,:) = onethree*q(0,:)+fivesixth*q(-1,:)-onesixth*q(-2, :)

    IS(1, :) = onefour*(q(2,:)-4.d0*q(1,:)+3.d0*q(0,:))**2 + &
               thirtwelve*(q(2,:)-2.d0*q(1,:)+q(0,:))**2

    IS(2, :) = onefour*(q(1,:)-q(-1,:))**2 + &
               thirtwelve*(q(1,:)-2.d0*q(0,:)+q(-1,:))**2
               
    IS(3, :) = onefour*(3.d0*q(0,:) - 4.d0*q(-1,:) +q(-2, :))**2 + &
               thirtwelve*(q(0,:)-2.d0*q(-1,:)+q(-2,:))**2

    alpha(1, :) = c(1)/(eps + IS(1,:))**p
    alpha(2, :) = c(2)/(eps + IS(2,:))**p
    alpha(3, :) = c(3)/(eps + IS(3,:))**p

    omega(1,:) = alpha(1, :)/(sum(alpha, dim=1))
    omega(2,:) = alpha(2, :)/(sum(alpha, dim=1))
    omega(3,:) = alpha(3, :)/(sum(alpha, dim=1))
    qR = omega(1,:)*flux(1,:) + omega(2,:)*flux(2,:) + omega(3,:)*flux(3,:)
  end subroutine weno_5th

  pure function minmod(a, b) result(res)
    real(rt) , intent(in) :: a, b
    real(rt) :: res
    if (a*b> 0) then
      if (abs(a)>abs(b)) then
        res = b
      else
        res = a
      end if
    else
      res = 0.d0
    endif
  end function minmod

end module advection_module
