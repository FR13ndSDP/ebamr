#ifndef _NC_F_H_
#define _NC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void initdata_f(const int* level, const amrex_real* time,
                    const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_3D(state),
                    const amrex_real* dx, const amrex_real* problo);


    void nc_tag_dengrad(const int* lo, const int* hi,
                           char* tag, const int* tlo, const int* thi,
                           amrex_real* rho, const int* rlo, const int* rhi,
                           const void* flag, const int* flo, const int* fhi,
                           const amrex_real* dengrad, const char* tagval, const char* clearval);

    void nullfill (BL_FORT_FAB_ARG_3D(state),
                       const int* dlo, const int* dhi,
                       const amrex_real* dx, const amrex_real* glo,
                       const amrex_real* time, const int* bc);

    void nc_denfill (BL_FORT_FAB_ARG_3D(state),
                    const int* dlo, const int* dhi,
                    const amrex_real* dx, const amrex_real* glo,
                    const amrex_real* time, const int* bc);
    void nc_hypfill (BL_FORT_FAB_ARG_3D(state),
                    const int* dlo, const int* dhi,
                    const amrex_real* dx, const amrex_real* glo,
                    const amrex_real* time, const int* bc);

    void nc_init_fort(const int* phys_bc_lo, const int* phys_bc_hi,
                      int interior, int inflow, int outflow, int symmetry, int slipwall, int noslipwall,
                      int myproc, const amrex_real* problo, const amrex_real* probhi);

    void nc_derpres (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                      const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const amrex_real* delta, const amrex_real* xlo,
                      const amrex_real* time, const amrex_real* dt, const int* bcrec,
                      const int* level, const int* grid_no);

    void nc_dervel (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex_real* delta, const amrex_real* xlo,
                     const amrex_real* time, const amrex_real* dt, const int* bcrec,
                     const int* level, const int* grid_no);
#ifdef __cplusplus
}
#endif

#endif