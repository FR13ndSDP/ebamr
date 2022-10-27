module nc_dudt_module

    use amrex_fort_module, only : rt=>amrex_real
    use nc_module, only : gamma, cv, smallp, smallr, nvar, qvar, qrho, qu, &
    qv, qw, qp, urho, ueden, umx, umy, umz, do_reflux, do_visc, refactor_scheme, NND, nghost_plm
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    implicit none

    private

    public:: compute_dudt, compute_divop, c2prim
  
  contains
  
    subroutine compute_dudt (lo,hi, dudt, utlo, uthi, &
         u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx,dt,level) &
         bind(c,name='compute_dudt')
      use advection_module, only : compute_flux
      use viscous_module, only : compute_viscous_flux
      integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi,fxlo,fxhi,fylo,fyhi,fzlo,fzhi
      real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
      real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
      real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),qvar)
      real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),qvar)
      real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),qvar)

      real(rt), intent(in) :: dx(3), dt
      integer, intent(in) :: level
  
      integer :: qlo(3), qhi(3)
      real(rt), dimension(:,:,:,:), pointer, contiguous :: q, fhx, fhy, fhz, fdx, fdy, fdz

      ! primitive is grown by nghost_plm
      qlo = lo - nghost_plm
      qhi = hi + nghost_plm
      call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
      call amrex_allocate(fhx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,qvar)
      call amrex_allocate(fhy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,qvar)
      call amrex_allocate(fhz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,qvar)
      call amrex_allocate(fdx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,qvar)
      call amrex_allocate(fdy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,qvar)
      call amrex_allocate(fdz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,qvar)

      call c2prim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
  
      ! compute face flux
      call compute_flux(q, qlo, qhi, lo, hi, dx, fhx, fhy, fhz)
      fx = fhx
      fy = fhy
      fz = fhz

      ! compute viscous flux
      if (do_visc) then
        call compute_viscous_flux(q, qlo, qhi, lo, hi, dx, fdx, fdy, fdz)
        fx = fx + fdx
        fy = fy + fdy
        fz = fz + fdz
        ! print *, "visc percent :", max(maxval(abs(fdx))/maxval(abs(fx)), &
        !                                max(maxval(abs(fdy))/maxval(abs(fy)), &
        !                                maxval(abs(fdz))/maxval(abs(fz))))
      end if

      ! get rhs by summing flux across face (div)
      call compute_divop (lo,hi,qvar,dx,dudt,utlo,uthi, &
           fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dt, level)
  
      call amrex_deallocate(q)
      call amrex_deallocate(fhx)
      call amrex_deallocate(fhy)
      call amrex_deallocate(fhz)
      call amrex_deallocate(fdx)
      call amrex_deallocate(fdy)
      call amrex_deallocate(fdz)
    end subroutine compute_dudt

    subroutine compute_divop(lo,hi,ncomp,dx,ut,utlo,uthi, &
        fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dt, level)
        integer, intent(in) :: lo(3),hi(3),ncomp,utlo(3),uthi(3),fxlo(3),fxhi(3), fylo(3),fyhi(3),fzlo(3),fzhi(3)
        real(rt), intent(inout) :: ut( utlo(1): uthi(1), utlo(2): uthi(2), utlo(3): uthi(3), ncomp)
        real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),ncomp)
        real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),ncomp)
        real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),ncomp)
        real(rt), intent(in) :: dx(3), dt
        integer, intent(in) :: level
    
        integer :: i,j,k,n
        real(rt) :: dxinv(3)
    
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
    end subroutine compute_divop

    subroutine c2prim(lo, hi, u, ulo, uhi, q, qlo, qhi)
      integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), qlo(3), qhi(3)
      real(rt), intent(in) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3), nvar)
      real(rt), intent(inout) :: q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3), qvar)
  
      integer :: i,j,k
      real(rt) :: rhoinv, gm1 = gamma - 1.d0
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            q(i,j,k,qrho) = max(u(i,j,k,urho), smallr)
            rhoinv = 1.d0/q(i,j,k,qrho)
            q(i,j,k,qu) = u(i,j,k,umx) * rhoinv
            q(i,j,k,qv) = u(i,j,k,umy) * rhoinv
            q(i,j,k,qw) = u(i,j,k,umz) * rhoinv
            q(i,j,k,qp) = gm1 * &
              (u(i,j,k,ueden) - 0.5d0*q(i,j,k,qrho)*(q(i,j,k,qu)**2 + q(i,j,k,qv)**2 + q(i,j,k,qw)**2))
          end do
        end do
      end do
    end subroutine c2prim

end module nc_dudt_module