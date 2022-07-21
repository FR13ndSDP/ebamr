module fluxsplit_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: flux_split

contains

  ! vanLeer
  subroutine flux_split(primL, primR, flux)
    use nc_module, only : gamma, qrho, qu, qv, qw, qp, qvar
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    real(rt), intent(in) :: primL(qvar), primR(qvar)
    real(rt), intent(inout) :: flux(qvar)
    real(rt) :: fp(qvar), fn(qvar)
    real(rt) :: cL, cR, ML, MR, Mp, tmp, Mn, tmpn
    real(rt) :: rhoL, rhoR, uL, uR, vL, vR, wL, wR, pL, pR

    real(rt) :: gm1 = gamma - 1.d0

    rhoL = primL(qrho); uL = primL(qu); vL = primL(qv); wL = primL(qw); pL = primL(qp);
    rhoR = primR(qrho); uR = primR(qu); vR = primR(qv); wR = primR(qw); pR = primR(qp);

    cL = sqrt(gamma * pL/rhoL)
    cR = sqrt(gamma * pR/rhoR)

    ML = uL/cL
    MR = uR/cR

    if (ML >= 1.d0) then
      fp(1) = rhoL * uL;
      fp(2) = rhoL * uL * uL + pL;
      fp(3) = rhoL * uL * vL;
      fp(4) = rhoL * uL * wL;
      fp(5) = uL * (gamma * pL / gm1 + 0.5d0 * rhoL * (uL * uL + vL * vL + wL*wL));
    else if (abs(ML) < 1.0) then
      Mp = 0.25d0 * (1.d0 + ML) * (1.d0 + ML);
      tmp = rhoL * cL * Mp;
      fp(1) = tmp;
      fp(2) = tmp * (gm1 * uL + 2.d0 * cL) / gamma;
      fp(3) = tmp * vL;
      fp(3) = tmp * wL;
      fp(5) = tmp * ((gm1 * uL + 2.d0 * cL)**2 * 0.5d0 / (gamma**2 - 1.d0) &
        + 0.5d0 * (vL * vL + wL * wL));
    else if (ML <= -1.d0) then
      fp = 0.d0
    end if

    if (MR >= 1.d0) then
      fn = 0.d0
    else if (abs(MR) < 1.d0) then
      Mn = -0.25d0 * (MR - 1.d0) * (MR - 1.d0);
      tmpn = rhoR * cR * Mn;
      fn(1) = tmpn;
      fn(2) = tmpn * (gm1 * uR - 2.d0 * cR) / gamma;
      fn(3) = tmpn * vR;
      fn(4) = tmpn * wR;
      fn(5) = tmpn * ((gm1 * uR - 2.d0 * cR)**2 * 0.5d0 / (gamma**2 - 1.d0) &
        + 0.5d0 * (vR * vR + wR * wR));
    else if (MR <= -1.d0) then
      fn(1) = rhoR * uR;
      fn(2) = rhoR * uR * uR + pR;
      fn(3) = rhoR * uR * vR;
      fn(4) = rhoR * uR * wR;
      fn(5) = uR * (gamma * pR / gm1 + 0.5d0 * rhoR * (uR * uR + vR * vR + wR*wR));
    end if
    flux = fp + fn

  end subroutine flux_split
end module fluxsplit_module
