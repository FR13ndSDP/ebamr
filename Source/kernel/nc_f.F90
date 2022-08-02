module nc_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  ! flux split scheme control
  integer, save, public :: split_scheme = 3
  integer, parameter, public :: vanLeer = 1
  integer, parameter, public :: HLLC = 2
  integer, parameter, public :: SW = 3
  integer, parameter, public :: Roe = 4
  
  ! L R flux scheme
  integer, save, public :: refactor_scheme = 1
  integer, parameter, public :: NND = 1
  integer, save, public :: nghost_plm = 2
  ! do reflux or not
  logical, save, public :: do_reflux = .false.
  
  ! do viscous term computing or not
  logical, save, public :: do_visc = .true.

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

  real(rt), parameter, public :: gamma = 1.4d0, cv = 717.5d0, R = 287.d0, Pr = 0.72d0
  
  public :: init_fort

contains

  subroutine init_fort() bind(c,name='init_fort')
    use amrex_parmparse_module
    integer :: i_do_reflux = 0, i_do_visc = 1
    type(amrex_parmparse) :: pp

    call amrex_parmparse_build(pp,"nc")

    call pp%query("split_scheme", split_scheme)
    call pp%query("refactor_scheme", refactor_scheme)
    call pp%query("do_reflux", i_do_reflux)
    call pp%query("do_visc", i_do_visc)
    do_reflux = i_do_reflux .ne. 0
    do_visc = i_do_visc .ne. 0

    if (refactor_scheme .ne. 1) nghost_plm = 3
    call amrex_parmparse_destroy(pp)

  end subroutine init_fort
end module nc_module
