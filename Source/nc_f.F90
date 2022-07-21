module nc_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  ! these flags must be the same as in CNS.H
  integer, parameter, public :: levmsk_interior   = 0 ! valid cells
  integer, parameter, public :: levmsk_covered    = 1 ! ghost cells covered by valid cells of this level
  integer, parameter, public :: levmsk_notcovered = 2 ! ghost cells not covered
  integer, parameter, public :: levmsk_physbnd    = 3 ! outside domain

  ! conservative variables index
  integer, parameter, public :: URHO  = 1
  integer, parameter, public :: UMX   = 2
  integer, parameter, public :: UMY   = 3
  integer, parameter, public :: UMZ   = 4
  integer, parameter, public :: UEDEN = 5
  integer, parameter, public :: NVAR  = 5

  ! primitive variables index
  integer, parameter, public :: QRHO   = 1
  integer, parameter, public :: QU     = 2
  integer, parameter, public :: QV     = 3
  integer, parameter, public :: QW     = 4
  integer, parameter, public :: QP     = 5
  integer, parameter, public :: QVAR   = 5


  real(rt), parameter, public :: smallp = 1.d-10
  real(rt), parameter, public :: smallr = 1.d-19

  ! boundary condition information
  integer, save, public :: physbc_lo(3)
  integer, save, public :: physbc_hi(3)
  integer, save, public :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! problem domain
  real(rt), save, public :: problo(3), probhi(3), center(3)

  integer, save, public :: myproc

  real(rt), parameter, public :: gamma = 1.4d0, cv = 717.5d0, R = 287.d0
  
  public :: nc_init_fort

contains

  subroutine nc_init_fort (physbc_lo_in, physbc_hi_in, &
       Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in, &
       myproc_in, problo_in, probhi_in) &
       bind(c,name='nc_init_fort')
    use amrex_parmparse_module
    use amrex_eb_flux_reg_nd_module, only : amrex_eb_disable_reredistribution
    integer, intent(in) :: physbc_lo_in(3), physbc_hi_in(3)
    integer, value, intent(in) :: Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
    integer, value, intent(in) :: myproc_in
    real(rt), intent(in) :: problo_in(3), probhi_in(3)

    physbc_lo = physbc_lo_in
    physbc_hi = physbc_hi_in

    Interior   = Interior_in
    Inflow     = Inflow_in
    Outflow    = Outflow_in
    Symmetry   = Symmetry_in
    SlipWall   = SlipWall_in
    NoSlipWall = NoSlipWall_in

    myproc = myproc_in

    problo = problo_in
    probhi = probhi_in
    center = 0.5d0*(problo+probhi)

  end subroutine nc_init_fort

end module nc_module
