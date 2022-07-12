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
#ifdef __cplusplus
}
#endif

#endif