module nc_ibm_module
    use amrex_fort_module, only : rt=>amrex_real
    use nc_module, only : gamma, cv, smallp, smallr, nvar, qvar, qrho, qu, &
    qv, qw, qp, urho, ueden, umx, umy, umz, nghost_plm
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate

    implicit none
    private

    public :: eb_compute_dudt, eb_compute_divop

contains

    subroutine eb_compute_divop(lo,hi,ncomp,dx,ut,utlo,uthi, &
        fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dt, level)
        integer, intent(in) :: lo(3),hi(3),ncomp,utlo(3),uthi(3),fxlo(3),fxhi(3), fylo(3),fyhi(3),fzlo(3),fzhi(3)
        real(rt), intent(inout) :: ut( utlo(1): uthi(1), utlo(2): uthi(2), utlo(3): uthi(3), ncomp)
        real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),ncomp)
        real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),ncomp)
        real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),ncomp)
        real(rt), intent(in) :: dx(3), dt
        integer, intent(in) :: level
    
        integer :: i,j,k,n
        real(rt) :: dxinv(3), coeff
    
        dxinv = 1.d0/dx
    
        do       k = lo(3),hi(3)
          do     j = lo(2),hi(2)
              do i = lo(1),hi(1)
                ut(i,j,k,:) = (fx(i,j,k,:)-fx(i+1,j,k,:)) * dxinv(1) &   
                      +       (fy(i,j,k,:)-fy(i,j+1,k,:)) * dxinv(2) &
                      +       (fz(i,j,k,:)-fz(i,j,k+1,:)) * dxinv(3)
              end do
          end do
        end do

    end subroutine eb_compute_divop

    subroutine eb_compute_dudt (lo,hi, dudt, utlo, uthi, &
      u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,flag,fglo,fghi, &
      volfrac,vlo,vhi, bcent,blo,bhi, &
      apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
      centx,cxlo,cxhi,centy,cylo,cyhi,centz,czlo,czhi, &
      as_crse_in, rr_drho_crse, rdclo, rdchi, rr_flag_crse, rfclo, rfchi, &
      as_fine_in, dm_as_fine, dflo, dfhi, &
      levmsk, lmlo, lmhi, &
      dx,dt,level) &
      bind(c,name='eb_compute_dudt')
   use nc_dudt_module, only : c2prim
   use advection_module, only : compute_flux
   integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi, &
        vlo,vhi,axlo,axhi,aylo,ayhi,azlo,azhi, &
        cxlo,cxhi,cylo,cyhi,czlo,czhi, &
        fglo,fghi, blo, bhi, fxlo,fxhi, fylo,fyhi, fzlo,fzhi, &
        rdclo, rdchi, rfclo, rfchi, dflo, dfhi, lmlo, lmhi
   integer, intent(in) :: as_crse_in, as_fine_in
   real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
   real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
   real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
   real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
   real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
   integer , intent(in) ::  flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
   real(rt), intent(in) :: volfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
   real(rt), intent(in) :: bcent  (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
   real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
   real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
   real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
   real(rt), intent(in) :: centx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
   real(rt), intent(in) :: centy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
   real(rt), intent(in) :: centz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
   real(rt), intent(inout) :: rr_drho_crse(rdclo(1):rdchi(1),rdclo(2):rdchi(2),rdclo(3):rdchi(3),nvar)
   integer,  intent(in) ::  rr_flag_crse(rfclo(1):rfchi(1),rfclo(2):rfchi(2),rfclo(3):rfchi(3))
   real(rt), intent(out) :: dm_as_fine(dflo(1):dfhi(1),dflo(2):dfhi(2),dflo(3):dfhi(3),nvar)
   integer,  intent(in) ::  levmsk (lmlo(1):lmhi(1),lmlo(2):lmhi(2),lmlo(3):lmhi(3))
   real(rt), intent(in) :: dx(3), dt
   integer, intent(in) :: level
 
     integer :: qlo(3), qhi(3)
     real(rt), dimension(:,:,:,:), pointer, contiguous :: q, fhx, fhy, fhz, fdx, fdy, fdz

     ! primitive is grown by nghost_plm
     qlo = lo - nghost_plm
     qhi = hi + nghost_plm
     call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
     call amrex_allocate(fdx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,qvar)
     call amrex_allocate(fdy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,qvar)
     call amrex_allocate(fdz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,qvar)

     call c2prim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
 
     ! compute face flux
     call compute_flux(q, qlo, qhi, lo, hi, dx, fx, fy, fz)

     ! get rhs by summing flux across face (div)
     call eb_compute_divop (lo,hi,qvar,dx,dudt,utlo,uthi, &
          fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dt, level)
 
     call amrex_deallocate(q)
     call amrex_deallocate(fdx)
     call amrex_deallocate(fdy)
     call amrex_deallocate(fdz)
    end subroutine eb_compute_dudt

end module nc_ibm_module