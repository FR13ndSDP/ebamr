
#include <NC.H>
#include <NC_F.H>

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_MultiCutFab.H>
#include <climits>

using namespace amrex;

#if __cplusplus < 201703L
constexpr int NC::level_mask_interior;
constexpr int NC::level_mask_covered;
constexpr int NC::level_mask_notcovered;
constexpr int NC::level_mask_physbnd;
constexpr int NC::NUM_GROW;
#endif

BCRec     NC::phys_bc;

int       NC::verbose = 0;
IntVect   NC::hydro_tile_size {AMREX_D_DECL(1024,16,16)};
Real      NC::cfl       = 0.3;
int       NC::do_reflux = 1;
int       NC::do_gravity = 0;
int       NC::refine_max_dengrad_lev   = -1;
Real      NC::refine_dengrad           = 1.0e10;
std::string NC::time_integration       = "RK2";
Vector<RealBox> NC::refine_boxes;

NC::NC ()
{}

NC::NC (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if (do_reflux && level > 0) {
        flux_reg.define(bl, papa.boxArray(level-1),
                        dm, papa.DistributionMap(level-1),
                        level_geom, papa.Geom(level-1),
                        papa.refRatio(level-1), level, NUM_STATE);
    }

    buildMetrics();
}

NC::~NC ()
{}

void
NC::init (AmrLevel& old)
{
    auto& oldlev = dynamic_cast<NC&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
}

void
NC::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

void
NC::initData ()
{
    BL_PROFILE("NC::initData()");

    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real* prob_hi = geom.ProbHi();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.validbox();
      initdata_f(&level, &cur_time,
                   BL_TO_FORTRAN_BOX(box),
                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                   dx, prob_lo, prob_hi);
    }
}

void NC::buildMetrics()
{
    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    if (std::abs(dx[0]-dx[1]) > 1.e-12*dx[0] || std::abs(dx[0]-dx[2]) > 1.e-12*dx[0]) {
        amrex::Abort("CNS: must have dx == dy == dz\n");
    }

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

    volfrac = &(ebfactory.getVolFrac());
    bndrycent = &(ebfactory.getBndryCent());
    areafrac = ebfactory.getAreaFrac();
    facecent = ebfactory.getFaceCent();

    level_mask.clear();
    level_mask.define(grids,dmap,1,1);
    level_mask.BuildMask(geom.Domain(), geom.periodicity(),
                         level_mask_covered,
                         level_mask_notcovered,
                         level_mask_physbnd,
                         level_mask_interior);
}

void
NC::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& ref_ratio,
                       Vector<Real>&          dt_level,
                       Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  if (level > 0) {
    return;
  }

  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
NC::computeNewDt (int                    finest_level,
                   int                    sub_cycle,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& ref_ratio,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) {
        return;
    }

    for (int i = 0; i <= finest_level; i++)
    {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1)
    {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }
    }
    else
    {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
NC::post_regrid (int lbase, int new_finest)
{
}

void
NC::post_timestep (int iteration)
{
    if (do_reflux && level < parent->finestLevel()) {
        NC& fine_level = getLevel(level+1);
        MultiFab& S_crse = get_new_data(State_Type);
        MultiFab& S_fine = fine_level.get_new_data(State_Type);
        fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac);
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }
}

void
NC::postCoarseTimeStep (Real time)
{
    // This only computes sum on level 0
    if (verbose >= 2) {
        printTotal();
    }
}

void
NC::printTotal () const
{
    const MultiFab& S_new = get_new_data(State_Type);
    MultiFab mf(grids, dmap, 1, 0);
    std::array<Real,5> tot;
    for (int comp = 0; comp < 5; ++comp) {
        MultiFab::Copy(mf, S_new, comp, 0, 1, 0);
        tot[comp] = mf.sum(0,true) * geom.ProbSize();
    }
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(tot.data(), 5, ParallelDescriptor::IOProcessorNumber());
            amrex::Print().SetPrecision(17) << "\n[NC] Total mass       is " << tot[0] << "\n"
                                            <<   "      Total x-momentum is " << tot[1] << "\n"
                                            <<   "      Total y-momentum is " << tot[2] << "\n"
                                            <<   "      Total z-momentum is " << tot[3] << "\n"
                                            <<   "      Total energy     is " << tot[4] << "\n";
#ifdef BL_LAZY
        });
#endif
}

void
NC::post_init (Real)
{
    if (level > 0) return;
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    if (verbose >= 2) {
        printTotal();
    }
}

void
NC::post_restart ()
{
}

void
NC::errorEst (TagBoxArray& tags, int, int, Real time, int, int)
{
    BL_PROFILE("NC::errorEst()");

    if (!refine_boxes.empty())
    {
        const Real* problo = geom.ProbLo();
        const Real* dx = geom.CellSize();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(tags); mfi.isValid(); ++mfi)
        {
            auto& fab = tags[mfi];
            const Box& bx = fab.box();
            for (BoxIterator bi(bx); bi.ok(); ++bi)
            {
                const IntVect& cell = bi();
                RealVect pos {AMREX_D_DECL((cell[0]+0.5)*dx[0]+problo[0],
                                           (cell[1]+0.5)*dx[1]+problo[1],
                                           (cell[2]+0.5)*dx[2]+problo[2])};
                for (const auto& rbx : refine_boxes) {
                    if (rbx.contains(pos)) {
                        fab(cell) = TagBox::SET;
                    }
                }
            }
        }
    }

    if (level < refine_max_dengrad_lev)
    {
        int ng = 1;
        const auto& rho = derive("rho", time, ng);
        const MultiFab& S_new = get_new_data(State_Type);

        const char   tagval = TagBox::SET;
        const char clearval = TagBox::CLEAR;

        auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S_new.Factory());
        auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rho,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& flag = flags[mfi];

            if (flag.getType(bx) != FabType::covered)
            {
                nc_tag_dengrad(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(tags[mfi]),
                                    BL_TO_FORTRAN_ANYD((*rho)[mfi]),
                                    BL_TO_FORTRAN_ANYD(flag),
                                    &refine_dengrad, &tagval, &clearval);
            }
        }
    }
}

void
NC::read_params ()
{
    ParmParse pp("nc");

    pp.query("v", verbose);

    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

    pp.query("cfl", cfl);

    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    pp.query("do_reflux", do_reflux);

    pp.query("do_gravity", do_gravity);

    pp.query("refine_max_dengrad_lev", refine_max_dengrad_lev);
    pp.query("refine_dengrad", refine_dengrad);

    int irefbox = 0;
    Vector<Real> refboxlo, refboxhi;
    while (pp.queryarr(("refine_box_lo_"+std::to_string(irefbox)).c_str(), refboxlo))
    {
        pp.getarr(("refine_box_hi_"+std::to_string(irefbox)).c_str(), refboxhi);
        refine_boxes.emplace_back(refboxlo.data(), refboxhi.data());
        ++irefbox;
    }
    pp.query("time_integration", time_integration);
}

void
NC::avgDown ()
{
    BL_PROFILE("NC::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level+1);

    MultiFab& S_crse =          get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
    volume.setVal(1.0);
    amrex::EB_average_down(S_fine, S_crse, volume, fine_lev.volFrac(),
                           0, S_fine.nComp(), fine_ratio);
}

Real NC::estTimeStep()
{
    BL_PROFILE("NC::estTimeStep()");

    Real estdt = std::numeric_limits<Real>::max();

    const Real* dx = geom.CellSize();
    const MultiFab& S = get_new_data(State_Type);

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:estdt)
#endif
    {
        Real dt = std::numeric_limits<Real>::max();
        for (MFIter mfi(S, true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            const auto& flag = flags[mfi];
            
            if (flag.getType(bx) != FabType::covered)
            {
                nc_estdt(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_ANYD(S[mfi]),
                            dx, &dt);
            }
        }
        estdt = std::min(estdt, dt);
    }

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);
    return estdt;
}

Real NC::initialTimeStep()
{
    return estTimeStep();
}

// core functions here
// S_new with 5 ghost cells
// Sborder with 5 ghost cells
// dSdt with 0 ghost cells
Real NC::advance(Real time, Real dt, int iteration, int /*ncycle*/)
{
    BL_PROFILE("NC::advance");

    for (int i=0; i<NUM_STATE_TYPE; i++)
    {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    // rhs
    MultiFab dSdt(grids, dmap, NUM_STATE, 0, MFInfo(), Factory());
    // state with ghost cell
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW, MFInfo(), Factory());

    EBFluxRegister* fine = nullptr; 
    EBFluxRegister* current = nullptr; 

    if (do_reflux && level < parent->finestLevel())
    {
        NC& fine_level = getLevel(level+1);
        fine = &fine_level.flux_reg;
        fine->reset();
    }

    if (do_reflux && level > 0)
    {
        current = &flux_reg;
    }

    // TODO: RK3 is not conservative in AMR
    // TODO: Use implicit time integration
    // TODO: use spectral deferred correction method
    if (time_integration == "RK2")
    {
        // RK2 stage 1
        FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
        compute_dSdt(Sborder, dSdt, 0.5*dt, fine, current);
        // U^* = U^n + dt * dUdt^n
        // S_new = 1 * Sborder + dt * dSdt
        MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);

        // Rk2 stage 2
        // after fillpatch Sborder is U^*
        FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, NUM_STATE);
        compute_dSdt(Sborder, dSdt, 0.5*dt, fine, current);
        // S_new = 0.5 * U^* + 0.5 * U^n + 0.5*dt*dUdt^*
        MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NUM_STATE, 0);
        MultiFab::Saxpy(S_new, 0.5*dt, dSdt, 0, 0, NUM_STATE, 0);
    } else
    {
        // RK3 stage 1
        FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
        compute_dSdt(Sborder, dSdt, Real(dt/6.0), fine, current);
        // U^* = U^n + dt * dUdt^n
        // S_new = 1 * Sborder + dt * dSdt
        MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);

        // Rk3 stage 2
        // after fillpatch Sborder is U^*
        FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, NUM_STATE);
        compute_dSdt(Sborder, dSdt, Real(dt/6.0), fine, current);
        // S_new = 0.25 * U^* + 0.75 * U^n + 0.25*dt*dUdt^*
        MultiFab::LinComb(S_new, 0.25, Sborder, 0, 0.75, S_old, 0, 0, NUM_STATE, 0);
        MultiFab::Saxpy(S_new, 0.25*dt, dSdt, 0, 0, NUM_STATE, 0);

        // Rk3 stage 3
        // after fillpatch Sborder is U^*
        FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, NUM_STATE);
        compute_dSdt(Sborder, dSdt, Real(2.0*dt/3.0), fine, current);
        MultiFab::LinComb(S_new, 2.0/3.0, Sborder, 0, 1.0/3.0, S_old, 0, 0, NUM_STATE, 0);
        MultiFab::Saxpy(S_new, 2.0/3.0*dt, dSdt, 0, 0, NUM_STATE, 0);
    }
    return dt;
}

// compute rhs flux
// bx with no ghost cell
// flux is defined on box face
void NC::compute_dSdt(const amrex::MultiFab &S, amrex::MultiFab &dSdt, amrex::Real dt, amrex::EBFluxRegister *fine, amrex::EBFluxRegister *current)
{
    BL_PROFILE("NC::compute_dSdt");

    const Real* dx = geom.CellSize();
    const int ncomp = dSdt.nComp();

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,AMREX_SPACEDIM> flux;
        FArrayBox dm_as_fine(Box::TheUnitBox(), ncomp);
        FArrayBox fab_drho_as_crse(Box::TheUnitBox(), ncomp);
        IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());

        for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true));
                    mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            const auto& flag = flags[mfi];

            if (flag.getType(bx) == FabType::covered) {
                dSdt[mfi].setVal<RunOn::Host>(0.0, bx , 0, ncomp);
            } else {
                // flux is used to store centroid flux needed for reflux
                // for flux_x in x direction is nodal, in other directions centroid
                // for flux_y in y ...
                for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
                    flux[idim].resize(amrex::surroundingNodes(bx,idim),ncomp);
                }

                // no cut cell around
                if (flag.getType(amrex::grow(bx,1)) == FabType::regular)
                {
                    compute_dudt(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                                BL_TO_FORTRAN_ANYD(S[mfi]),
                                BL_TO_FORTRAN_ANYD(flux[0]),
                                BL_TO_FORTRAN_ANYD(flux[1]),
                                BL_TO_FORTRAN_ANYD(flux[2]),
                                dx, &dt, &level);
                    
                    if (do_gravity != 0) {
                        const Real g = -9.8;
                        const int irho = Density;
                        const int imz = Zmom;
                        const int irhoE = Eden;
                        auto const& dsdtfab = dSdt.array(mfi); 
                        auto const& sfab = S.array(mfi); 
                        amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            dsdtfab(i,j,k,imz) += g * sfab(i,j,k,irho);
                            dsdtfab(i,j,k,irhoE) += g * sfab(i,j,k,imz);
                        });
                    }

                    // TODO: reflux for EB is too complicated!
                    if (do_reflux)
                    {
                        // the flux registers from the coarse or fine grid perspective
                        // NOTE: the flux register associated with flux_reg[lev] is associated
                        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
                        if (current) {
                        // update the lev/lev-1 flux register (index lev)
                            for (int i=0; i<AMREX_SPACEDIM; i++)
                                current->FineAdd(mfi, {&flux[0], &flux[1], &flux[2]}, dx, dt, RunOn::Cpu);
                        }

                        if (fine) {
                        // update the lev+1/lev flux register (index lev+1)
                            for (int i=0; i<AMREX_SPACEDIM; i++)
                                fine->CrseAdd(mfi, {&flux[0], &flux[1], &flux[2]}, dx, dt, RunOn::Cpu);
                        }
                    }
                }
                else
                {
                    // cut cells and its neighbor

                    FArrayBox* p_drho_as_crse = (fine) ?
                            fine->getCrseData(mfi) : &fab_drho_as_crse;
                    const IArrayBox* p_rrflag_as_crse = (fine) ?
                        fine->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                    if (current) {
                        dm_as_fine.resize(amrex::grow(bx,1),ncomp);
                    }

                    int as_fine = (fine != nullptr);
                    int as_crse = (current != nullptr);

                    (*areafrac[0])[mfi];
                    eb_compute_dudt(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                                    BL_TO_FORTRAN_ANYD(S[mfi]),
                                    BL_TO_FORTRAN_ANYD(flux[0]),
                                    BL_TO_FORTRAN_ANYD(flux[1]),
                                    BL_TO_FORTRAN_ANYD(flux[2]),
                                    BL_TO_FORTRAN_ANYD(flag),
                                    BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                                    &as_fine,
                                    BL_TO_FORTRAN_ANYD(*p_drho_as_crse),
                                    BL_TO_FORTRAN_ANYD(*p_rrflag_as_crse),
                                    &as_crse,
                                    BL_TO_FORTRAN_ANYD(dm_as_fine),
                                    BL_TO_FORTRAN_ANYD(level_mask[mfi]),
                                    dx, &dt,&level);

                    if (fine) {
                        fine->CrseAdd(mfi, {&flux[0],&flux[1],&flux[2]}, dx,dt,
                                            (*volfrac)[mfi],
                                            {&((*areafrac[0])[mfi]),
                                             &((*areafrac[1])[mfi]),
                                             &((*areafrac[2])[mfi])},
                                            RunOn::Cpu);
                    }

                    if (current) {
                        current->FineAdd(mfi, {&flux[0],&flux[1],&flux[2]}, dx,dt,
                                            (*volfrac)[mfi],
                                            {&((*areafrac[0])[mfi]),
                                             &((*areafrac[1])[mfi]),
                                             &((*areafrac[2])[mfi])},
                                            dm_as_fine,
                                            RunOn::Cpu);
                    }
                }
            }
        }
    }
}