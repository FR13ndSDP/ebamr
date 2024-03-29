module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  real(rt), save :: p_l   = 1.d0
  real(rt), save :: p_r   = 1.d0
  real(rt), save :: rho_l = 1.d0
  real(rt), save :: rho_r = 1.d0
  real(rt), save :: u_l   = 3.d0*sqrt(1.4d0)
  real(rt), save :: v_l   = 0.d0
  real(rt), save :: u_r   = 3.d0*sqrt(1.4d0)
  real(rt), save :: v_r   = 0.d0
end module probdata_module


subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  use amrex_fort_module, only : rt => amrex_real
  use amrex_parmparse_module
  use probdata_module
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(*), probhi(*)
  type(amrex_parmparse) :: pp
  call amrex_parmparse_build(pp,"prob")
  call pp%query("p_l",p_l)
  call pp%query("p_r",p_r)
  call pp%query("rho_l",rho_l)
  call pp%query("rho_r",rho_r)
  call pp%query("u_l",u_l)
  call pp%query("u_r",u_r)
  call amrex_parmparse_destroy(pp)
end subroutine amrex_probinit


subroutine initdata_f(level, time, lo, hi, u, ulo, uhi, dx, prob_lo, prob_hi) bind(C, name="initdata_f")
  use amrex_fort_module, only : rt => amrex_real
  use nc_module, only : nvar, urho, umx, umy, umz, ueden, gamma, cv
  use probdata_module
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3), prob_hi(3)

  integer :: i,j,k
  real(rt) :: x, y, z, Pt, rhot, uxt, vxt

  do k = lo(3), hi(3)
    z = prob_lo(3) + (k+0.5d0)*dx(3)
    do j = lo(2), hi(2)
      y = prob_lo(2) + (j+0.5d0)*dx(2)
      do i = lo(1), hi(1)
        x = prob_lo(1) + (i+0.5d0)*dx(1)

        ! if ((x-0.5d0)**2+(y-0.5d0)**2+(z-0.5d0)**2 <= 0.04d0) then
        ! if (x<= 0.1d0) then
        !   Pt = p_l
        !   rhot = rho_l
        !   uxt = u_l
        !   vxt = v_l
        ! else
          Pt = p_l
          rhot = rho_l
          uxt = u_l
          vxt = v_l
        ! end if

        u(i,j,k,urho) = rhot
        u(i,j,k,umx) = rhot * uxt
        u(i,j,k,umy) = rhot * vxt
        u(i,j,k,umz) = 0.d0
        u(i,j,k,ueden) = Pt/(gamma - 1.d0) + 0.5d0*rhot*(uxt*uxt + vxt*vxt)
      end do
    end do
  end do

end subroutine initdata_f
