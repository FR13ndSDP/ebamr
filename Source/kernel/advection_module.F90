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
          else if (refactor_scheme .eq. 2) then
            call weno_5th(qL, qR, q(i-3:i+2,j,k,:))
          else
            call teno_5th(qL, qR, q(i-3:i+2,j,k,:))
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
          else if (refactor_scheme .eq. 2) then
            call weno_5th(qL, qR, q(i,j-3:j+2,k,:))
          else
            call teno_5th(qL, qR, q(i,j-3:j+2,k,:))
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
          else if (refactor_scheme .eq. 2) then
            call weno_5th(qL, qR, q(i,j,k-3:k+2,:))
          else
            call teno_5th(qL, qR, q(i,j,k-3:k+2,:))
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

  subroutine teno_5th(qL, qR, q)
    use nc_module, only : QVAR
    integer :: i
    real(rt), intent(inout) :: qL(qvar), qR(qvar)
    real(rt), intent(in) :: q(-3:2, qvar)
    real(rt) :: s11(qvar), s22(qvar), s33(qvar), s55(qvar), sumation(qvar), a1(qvar), &
                a2(qvar), a3(qvar), b1(qvar), b2(qvar), b3(qvar), v1(qvar), v2(qvar), &
                v3(qvar), w1(qvar), w2(qvar), w3(qvar)
    real(rt), parameter :: eps = 1d-40
    real(rt), parameter :: onesix = 1.d0/6.d0

    ! qL
    s11 = 13.d0*(q(-3,:)-2.d0*q(-2,:)+q(-1,:))**2 + 3.d0*(q(-3,:)-4.d0*q(-2,:)+3.d0*q(-1,:))**2
    s22 = 13.d0*(q(-2,:)-2.d0*q(-1,:)+q(0,:))**2 + 3.d0*(q(-2,:)-q(0,:))**2
    s33 = 13.d0*(q(-1,:)-2.d0*q(0,:)+q(1,:))**2 + 3.d0*(3.d0*q(-1,:)-4.d0*q(0,:)+q(1,:))**2

    s55 = abs(s11-s33)

    a1 = (1.d0+s55/(s11+eps))**6
    a2 = (1.d0+s55/(s22+eps))**6
    a3 = (1.d0+s55/(s33+eps))**6
    
    sumation = a1+a2+a3
    b1 = a1/sumation
    b2 = a2/sumation
    b3 = a3/sumation

    do i = 1,QVAR
      if (b1(i) < 1d-5) then
        b1(i) = 0.d0
      else
        b1(i) = 1.d0
      end if

      if (b2(i) < 1d-5) then
        b2(i) = 0.d0
      else
        b2(i) = 1.d0
      end if

      if (b3(i) < 1d-5) then
        b3(i) = 0.d0
      else
        b3(i) = 1.d0
      end if
    end do

    v1 = onesix*(2.d0*q(-3,:)-7.d0*q(-2,:)+5.d0*q(-1,:))
    v2 = onesix*(-q(-2,:)-q(-1,:)+2.d0*q(0,:))
    v3 = onesix*(-4.d0*q(-1,:)+5.d0*q(0,:)-q(1,:))

    a1 = 0.1d0*b1
    a2 = 0.6d0*b2
    a3 = 0.3d0*b3

    sumation = a1+a2+a3
    w1 = a1/sumation
    w2 = a2/sumation
    w3 = a3/sumation

    qL = q(-1,:)+w1*v1+w2*v2+w3*v3

    ! qR
    s11 = 13.d0*(q(2,:)-2.d0*q(1,:)+q(0,:))**2 + 3.d0*(q(2,:)-4.d0*q(1,:)+3.d0*q(0,:))**2
    s22 = 13.d0*(q(-1,:)-2.d0*q(0,:)+q(1,:))**2 + 3.d0*(q(1,:)-q(-1,:))**2
    s33 = 13.d0*(q(0,:)-2.d0*q(-1,:)+q(-2,:))**2 + 3.d0*(3.d0*q(0,:)-4.d0*q(-1,:)+q(-2,:))**2

    s55 = abs(s11-s33)

    a1 = (1.d0+s55/(s11+eps))**6
    a2 = (1.d0+s55/(s22+eps))**6
    a3 = (1.d0+s55/(s33+eps))**6
    
    sumation = a1+a2+a3
    b1 = a1/sumation
    b2 = a2/sumation
    b3 = a3/sumation

    do i = 1,QVAR
      if (b1(i) < 1d-5) then
        b1(i) = 0.d0
      else
        b1(i) = 1.d0
      end if

      if (b2(i) < 1d-5) then
        b2(i) = 0.d0
      else
        b2(i) = 1.d0
      end if

      if (b3(i) < 1d-5) then
        b3(i) = 0.d0
      else
        b3(i) = 1.d0
      end if
    end do

    v1 = onesix*(2.d0*q(2,:)-7.d0*q(1,:)+5.d0*q(0,:))
    v2 = onesix*(-q(1,:)-q(0,:)+2.d0*q(-1,:))
    v3 = onesix*(-4.d0*q(0,:)+5.d0*q(-1,:)-q(-2,:))

    a1 = 0.1d0*b1
    a2 = 0.6d0*b2
    a3 = 0.3d0*b3

    sumation = a1+a2+a3
    w1 = a1/sumation
    w2 = a2/sumation
    w3 = a3/sumation

    qR = q(0,:)+w1*v1+w2*v2+w3*v3
  end subroutine teno_5th

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
