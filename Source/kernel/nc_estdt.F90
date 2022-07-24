module nc_estdt_module
  use amrex_fort_module, only : rt=>amrex_real
  use nc_module, only : NVAR, gamma, urho, umx, umy, umz, ueden
  public :: nc_estdt
contains

  subroutine nc_estdt(lo, hi, u, ulo, uhi, dx, dt) bind(c, name='nc_estdt')
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3), nvar)
    real(rt), intent(inout) :: dt

    integer :: i,j,k
    real(rt) :: rhoinv, vx, vy, vz, p, cs

    do     k = lo(3),hi(3)
      do   j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhoinv = 1.d0/u(i,j,k,urho)
          vx = u(i,j,k,umx)*rhoinv
          vy = u(i,j,k,umy)*rhoinv
          vz = u(i,j,k,umz)*rhoinv
          p = (gamma-1.d0)*(u(i,j,k,ueden) - 0.5d0*u(i,j,k,urho)*(vx**2+vy**2+vz**2))
          cs = sqrt(gamma*p*rhoinv)
          dt = min(dt,dx(1)/(abs(vx)+cs),dx(2)/(abs(vy)+cs),dx(3)/(abs(vz)+cs))
        end do
      end do
    end do
  end subroutine nc_estdt

end module nc_estdt_module
