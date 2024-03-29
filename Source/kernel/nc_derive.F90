module nc_derive_module

  use amrex_fort_module, only : rt=>amrex_real
  use iso_fortran_env, only : stderr=>error_unit
  use nc_module, only : gamma, cv, smallp, smallr, nvar, qvar, qrho, qu, &
    qv, qw, qp, urho, ueden, umx, umy, umz
  implicit none
  private

  public :: nc_derpres, nc_dervel, nc_dertemp

contains

  ! This routine will derive pressure from internal energy density
subroutine nc_derpres(vel,v_lo,v_hi,nv, &
  dat,d_lo,d_hi,nc,lo,hi,domlo, &
  domhi,delta,xlo,time,dt,bc,level,grid_no) &
  bind(C, name="nc_derpres")
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: v_lo(3), v_hi(3), nv
  integer, intent(in) :: d_lo(3), d_hi(3), nc
  integer, intent(in) :: domlo(3), domhi(3)
  integer, intent(in) :: bc(3,2,nc)
  real(rt), intent(in) :: delta(3), xlo(3), time, dt
  real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  integer, intent(in) :: level, grid_no

  integer :: i, j, k
  real(rt) :: rhoinv, gm1 = gamma - 1.d0

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rhoinv = 1.d0/dat(i,j,k,1)
        vel(i,j,k) = (dat(i,j,k,5) - 0.5d0 * rhoinv * (dat(i,j,k,2)**2 + dat(i,j,k,3)**2+dat(i,j,k,4)**2))&
        * gm1
      end do
    end do
  end do

end subroutine nc_derpres


  ! This routine will derive the velocity from the momentum.
  subroutine nc_dervel(vel,v_lo,v_hi,nv, &
    dat,d_lo,d_hi,nc,lo,hi,domlo, &
    domhi,delta,xlo,time,dt,bc,level,grid_no) &
    bind(C, name="nc_dervel")
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3), nv
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    integer :: i, j, k

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          vel(i,j,k) = dat(i,j,k,2) / dat(i,j,k,1)
        end do
      end do
    end do

  end subroutine nc_dervel

  subroutine nc_dertemp(vel,v_lo,v_hi,nv, &
    dat,d_lo,d_hi,nc,lo,hi,domlo, &
    domhi,delta,xlo,time,dt,bc,level,grid_no) &
    bind(C, name="nc_dertemp")
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3), nv
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    integer :: i, j, k
    real(rt) :: rhoinv

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          rhoinv = 1.d0/dat(i,j,k,1)
          vel(i,j,k) = (dat(i,j,k,5) - 0.5d0 * rhoinv * (dat(i,j,k,2)**2 + dat(i,j,k,3)**2+dat(i,j,k,4)**2))&
          /(cv * dat(i,j,k,1))
        end do
      end do
    end do

  end subroutine nc_dertemp
end module nc_derive_module
