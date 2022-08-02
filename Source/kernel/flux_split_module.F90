module fluxsplit_module

  use amrex_fort_module, only : rt=>amrex_real
  use nc_module, only : split_scheme, vanLeer, HLLC, SW, Roe, gamma, qrho, qu, qv, qw, qp, qvar
  implicit none

  private :: flux_split_vanLeer, flux_split_HLLC, flux_split_SW, flux_split_Roe

  public :: flux_split
contains

  subroutine flux_split(primL, primR, flux)
    use amrex_error_module, only : amrex_error
    real(rt), intent(in) :: primL(qvar), primR(qvar)
    real(rt), intent(inout) :: flux(qvar)

    select case (split_scheme)
     case (vanLeer)
      call flux_split_vanLeer(primL, primR, flux)
     case (HLLC)
      call flux_split_HLLC(primL, primR, flux)
     case (SW)
      call flux_split_SW(primL, primR, flux)
     case (Roe)
      call flux_split_Roe(primL, primR, flux)
     case default
      call amrex_error("No other flux split schemes!")
    end select
  end subroutine flux_split

  ! vanLeer
  subroutine flux_split_vanLeer(primL, primR, flux)
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

    fp = 0.d0
    fn = 0.d0

    if (ML >= 1.d0) then
      fp(1) = rhoL * uL
      fp(2) = rhoL * uL * uL + pL
      fp(3) = rhoL * uL * vL
      fp(4) = rhoL * uL * wL
      fp(5) = uL * (gamma * pL / gm1 + 0.5d0 * rhoL * (uL * uL + vL * vL + wL*wL))
    else if (abs(ML) < 1.0) then
      Mp = 0.25d0 * (1.d0 + ML) * (1.d0 + ML)
      tmp = rhoL * cL * Mp
      fp(1) = tmp
      fp(2) = tmp * (gm1 * uL + 2.d0 * cL) / gamma
      fp(3) = tmp * vL
      fp(4) = tmp * wL
      fp(5) = tmp * ((gm1 * uL + 2.d0 * cL)**2 * 0.5d0 / (gamma**2 - 1.d0) &
        + 0.5d0 * (vL * vL + wL * wL))
    end if

    if (abs(MR) < 1.d0) then
      Mn = -0.25d0 * (MR - 1.d0) * (MR - 1.d0)
      tmpn = rhoR * cR * Mn
      fn(1) = tmpn
      fn(2) = tmpn * (gm1 * uR - 2.d0 * cR) / gamma
      fn(3) = tmpn * vR
      fn(4) = tmpn * wR
      fn(5) = tmpn * ((gm1 * uR - 2.d0 * cR)**2 * 0.5d0 / (gamma**2 - 1.d0) &
        + 0.5d0 * (vR * vR + wR * wR))
    else if (MR <= -1.d0) then
      fn(1) = rhoR * uR
      fn(2) = rhoR * uR * uR + pR
      fn(3) = rhoR * uR * vR
      fn(4) = rhoR * uR * wR
      fn(5) = uR * (gamma * pR / gm1 + 0.5d0 * rhoR * (uR * uR + vR * vR + wR*wR))
    end if
    flux = fp + fn

  end subroutine flux_split_vanLeer

  ! HLLC
  subroutine flux_split_HLLC(QL,QR,Flux)
    implicit none
    real(rt):: QL(5),QR(5),UL(5),UR(5),Flux(5)
    real(rt):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar    ! uu,vv,ww velocity
    real(rt):: p_pvrs,p_star,qql,qqr,Sl,Sr,Fl(5),Fr(5),S_star,tmpl,tmpr,tmp,F_starl(5),F_starr(5)

    dl=QL(1); uul=QL(2); vvl=QL(3);wwl=QL(4); pl=QL(5)
    dr=QR(1); uur=QR(2); vvr=QR(3);wwr=QR(4); pr=QR(5)
    UL(1)=dl; UL(2)=dl*uul; UL(3)=dl*vvl; UL(4)=dl*wwl; UL(5)=pl/(gamma-1.d0)+dl*(uul*uul+vvl*vvl+wwl*wwl)*0.5d0
    UR(1)=dr; UR(2)=dr*uur; UR(3)=dr*vvr; UR(4)=dr*wwr; UR(5)=pr/(gamma-1.d0)+dr*(uur*uur+vvr*vvr+wwr*wwr)*0.5d0
    al=sqrt(gamma*pl/dl); ar=sqrt(gamma*pr/dr)
    p_pvrs=0.5d0*(pl+pr)-(uur-uul)*(dl+dr)*(al+ar)*0.125d0
    p_star=max(0.d0,p_pvrs)

    if(p_star .le. pl) then
      qql=1
    else
      qql=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pl-1.d0) )
    endif
    if(p_star .le. pr) then
      qqr=1
    else
      qqr=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pr-1.d0) )
    endif
    Sl=uul-al*qql; Sr=uur+ar*qqr    ! speed of the lift and the right shockwaves

    Fl(1)=UL(2); Fl(2)=UL(2)*uul+pl; Fl(3)=UL(3)*uul; Fl(4)=UL(4)*uul; Fl(5)=uul*(UL(5)+pl)
    Fr(1)=UR(2); Fr(2)=UR(2)*uur+pr; Fr(3)=UR(3)*uur; Fr(4)=UR(4)*uur; Fr(5)=uur*(UR(5)+pr)
!----HLL---------------------------------

    if( Sl .ge. 0 ) then
      Flux=Fl
    else if (Sr .le. 0) then
      Flux=Fr
    else
      Flux=(Sr*Fl-Sl*Fr+Sl*Sr*(UR-UL))/(Sr-Sl)
    endif

  end subroutine flux_split_HLLC

  ! Steger-Warming
  subroutine flux_split_SW(QL,QR,Flux)
    implicit none
    real(rt):: QL(5),QR(5),Flux(5)
    real(rt):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar ,pr1   ! uu velocity
    real(rt):: tmp0,tmp1,tmp2,tmp3,E1P,E2P,E3P,E1M,E2M,E3M,fp(5),fm(5)


    dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5)
    dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
    al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
    ar=sqrt(gamma*pr/dr)  ! find a bug, removed

    tmp1=2.d0*(gamma-1.d0)
    tmp3=(3.d0-gamma)/(2.d0*(gamma-1.d0))
! eigenvalues---------
    E1P=(uul+abs(uul))*0.5d0
    E2P=(uul-al+abs(uul-al))*0.5d0
    E3P=(uul+al+abs(uul+al))*0.5d0
    tmp0=dl/(2.d0*gamma)
    fp(1)=tmp0*(tmp1*E1P+E2P+E3P)
    fp(2)=tmp0*(tmp1*E1P*uul+E2P*(uul-al)+E3P*(uul+al))
    fp(3)=tmp0*(tmp1*E1P*vvl+E2P*vvl+E3P*vvl)
    fp(4)=tmp0*(tmp1*E1P*wwl+E2P*wwl+E3P*wwl)
    fp(5)=tmp0*(E1P*(gamma-1.d0)*(uul*uul+vvl*vvl+wwl*wwl)+E2P*((uul-al)**2+vvl*vvl+wwl*wwl)*0.5d0 &
      +E3P*((uul+al)**2+vvl*vvl+wwl*wwl)*0.5d0+tmp3*al*al*(E2P+E3P))

    E1M=(uur-abs(uur))*0.5d0
    E2M=(uur-ar-abs(uur-ar))*0.5d0
    E3M=(uur+ar-abs(uur+ar))*0.5d0
    tmp0=dr/(2.d0*gamma)
    fm(1)=tmp0*(tmp1*E1M+E2M+E3M)
    fm(2)=tmp0*(tmp1*E1M*uur+E2M*(uur-ar)+E3M*(uur+ar))
    fm(3)=tmp0*(tmp1*E1M*vvr+E2M*vvr+E3M*vvr)
    fm(4)=tmp0*(tmp1*E1M*wwr+E2M*wwr+E3M*wwr)
    fm(5)=tmp0*(E1M*(gamma-1.d0)*(uur*uur+vvr*vvr+wwr*wwr)+E2M*((uur-ar)**2+vvr*vvr+wwr*wwr)*0.5d0 &
      +E3M*((uur+ar)**2+vvr*vvr+wwr*wwr)*0.5d0+tmp3*ar*ar*(E2M+E3M))
    Flux=fp+fm
  end subroutine flux_split_SW

  subroutine flux_split_Roe(QL,QR,Flux)
    implicit none
    real(rt):: QL(5),QR(5),Flux(5)
    real(rt):: dl,uul,vvl,wwl,pl,al, dr,uur,vvr,wwr,pr,ar    ! uu velocity
    real(rt):: avd,avu,avv,avw,avh,ava
    real(rt):: lamda(3),UL(5),UR(5),Fl(5),Fr(5)
    real(rt):: arr(5)
    real(rt):: hl,hr,tmpp,delt=0.1d0
    dl=QL(1); uul=QL(2); vvl=QL(3); wwl=QL(4); pl=QL(5)
    dr=QR(1); uur=QR(2); vvr=QR(3); wwr=QR(4); pr=QR(5)
    al=sqrt(gamma*pl/dl); ar=sqrt(gamma*pr/dr)
    hl=gamma*pl/((gamma-1.d0)*dl)+0.5*(uul*uul+vvl*vvl+wwl*wwl)
    hr=gamma*pr/((gamma-1.d0)*dr)+0.5*(uur*uur+vvr*vvr+wwr*wwr)

    UL(1)=dl; UL(2)=dl*uul; UL(3)=dl*vvl; UL(4)=dl*wwl; UL(5)=pl/(gamma-1.d0)+dl*(uul*uul+vvl*vvl+wwl*wwl)*0.5d0
    UR(1)=dr; UR(2)=dr*uur; UR(3)=dr*vvr; UR(4)=dr*wwr; UR(5)=pr/(gamma-1.d0)+dr*(uur*uur+vvr*vvr+wwr*wwr)*0.5d0
    Fl(1)=UL(2); Fl(2)=UL(2)*uul+pl; Fl(3)=UL(3)*uul; Fl(4)=UL(4)*uul; Fl(5)=uul*(UL(5)+pl)
    Fr(1)=UR(2); Fr(2)=UR(2)*uur+pr; Fr(3)=UR(3)*uur; Fr(4)=UR(4)*uur; Fr(5)=uur*(UR(5)+pr)

    tmpp=sqrt(dr/dl)
    avd=0.25d0*dl*(1.d0+tmpp)*(1.d0+tmpp)
    avu=(uul+tmpp*uur)/(1.d0+tmpp)
    avv=(vvl+tmpp*vvr)/(1.d0+tmpp)
    avw=(wwl+tmpp*wwr)/(1.d0+tmpp)
    avh=(hl+tmpp*hr)/(1.d0+tmpp)
    ava=sqrt((gamma-1.)*(avh-0.5d0*(avu*avu+avv*avv+avw*avw)))

    lamda(1)=abs(avu-ava)
    lamda(2)=abs(avu)
    lamda(3)=abs(avu+ava)

    if(lamda(1) < delt) then
      lamda(1)=( lamda(1)**2 +delt*delt)/(2.d0*delt)
    end if
    if(lamda(2) < delt) then
      lamda(2)=( lamda(2)**2 +delt*delt)/(2.d0*delt)
    end if
    if(lamda(3) < delt) then
      lamda(3)=( lamda(3)**2 +delt*delt)/(2.d0*delt)
    end if


    arr(2)=(UR(3)-UL(3))/(avv+1.e-8)-(dr-dl)
    arr(3)=(UR(4)-UL(4))/(avw+1.e-8)-(dr-dl)
    arr(4)=(gamma-1.d0)*((dr-dl)*(avh-avu*avu-avv*avv-avw*avw)+avu*(UR(2)-UL(2))+avv*(UR(3)-UL(3)) &
      +avw*(UR(4)-UL(4))-(UR(5)-UL(5)))/(ava*ava)
    arr(1)=((avu+ava)*(dr-dl)-UR(2)+UL(2)-ava*arr(4))/(2.d0*ava)
    arr(5)=dr-dl-(arr(1)+arr(4))
    flux(1)=0.5d0*(Fl(1)+Fr(1))-0.5d0*(lamda(1)*arr(1)+lamda(2)*arr(4)+lamda(3)*arr(5))
    flux(2)=0.5d0*(Fl(2)+Fr(2))-0.5d0*(lamda(1)*arr(1)*(avu-ava)+lamda(2)*arr(4)*avu+lamda(3)*arr(5)*(avu+ava))
    flux(3)=0.5d0*(Fl(3)+Fr(3))-0.5d0*(lamda(1)*arr(1)*avv+lamda(2)*arr(2)*avv+lamda(2)*arr(4)*avv+lamda(3)*arr(5)*avv)
    flux(4)=0.5d0*(Fl(4)+Fr(4))-0.5d0*(lamda(1)*arr(1)*avw+lamda(2)*arr(3)*avw+lamda(2)*arr(4)*avw+lamda(3)*arr(5)*avw)
    flux(5)=0.5d0*(Fl(5)+Fr(5))-0.5d0*(lamda(1)*arr(1)*(avh-avu*ava)+lamda(2)*arr(2)*avv*avv+lamda(2)*arr(3)*avw*avw  &
      +lamda(2)*arr(4)*0.5d0*(avu*avu+avv*avv+avw*avw)+lamda(3)*arr(5)*(avh+avu*ava))

  end subroutine flux_split_Roe
end module fluxsplit_module
