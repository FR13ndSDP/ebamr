module nc_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

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
  
end module nc_module
