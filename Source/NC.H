#ifndef _NC_H_
#define _NC_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBCellFlag.H>

using namespace amrex;
class NC : public amrex::AmrLevel
{
public:
    NC();
    NC(amrex::Amr &papa,
       int lev,
       const amrex::Geometry &level_geom,
       const amrex::BoxArray &bl,
       const amrex::DistributionMapping &dm,
       amrex::Real time);
    virtual ~NC();

    // disable copy
    NC(const NC &rhs) = delete;
    NC &operator=(const NC &rhs) = delete;

    // Restart from a checkpoint file.
    virtual void restart(amrex::Amr &papa,
                         std::istream &is,
                         bool bReadSpecial = false) override;

    // Write checkpoint
    virtual void checkPoint(const std::string &dir,
                            std::ostream &os,
                            amrex::VisMF::How how = amrex::VisMF::NFiles,
                            bool dump_old = true) override;

    // virtual std::string thePlotFileType () const override {
    //     return {"HyperCLaw-V1.1"};
    // }

    // Write state for debugging
    void printState(const MultiFab& mf);

    // Write a plotfile to specified directory.
    virtual void writePlotFile(const std::string &dir,
                               std::ostream &os,
                               amrex::VisMF::How how) override;

    // Initialize data on this level from another NC (during regrid).
    virtual void init(amrex::AmrLevel &old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init() override;

    // Initialize grid data at problem start-up.
    virtual void initData() override;

    // Advance grids at this level in time.
    virtual amrex::Real advance(amrex::Real time,
                                amrex::Real dt,
                                int iteration,
                                int ncycle) override;

    virtual void computeInitialDt(int finest_level,
                                  int sub_cycle,
                                  amrex::Vector<int> &n_cycle,
                                  const amrex::Vector<amrex::IntVect> &ref_ratio,
                                  amrex::Vector<amrex::Real> &dt_level,
                                  amrex::Real stop_time) override;

    virtual void computeNewDt(int finest_level,
                              int sub_cycle,
                              amrex::Vector<int> &n_cycle,
                              const amrex::Vector<amrex::IntVect> &ref_ratio,
                              amrex::Vector<amrex::Real> &dt_min,
                              amrex::Vector<amrex::Real> &dt_level,
                              amrex::Real stop_time,
                              int post_regrid_flag) override;

    virtual void post_regrid(int lbase, int new_finest) override;

    // Do work after timestep().
    virtual void post_timestep(int iteration) override;

    // After a full time step
    virtual void postCoarseTimeStep(amrex::Real time) override;

    // Do work after init().
    virtual void post_init(amrex::Real stop_time) override;

    virtual void post_restart() override;

    // Error estimation for regridding.
    virtual void errorEst(amrex::TagBoxArray &tb,
                          int clearval,
                          int tagval,
                          amrex::Real time,
                          int n_error_buf = 0,
                          int ngrow = 0) override;

    static int numGrow() { return NUM_GROW; }
    
    const amrex::MultiFab& volFrac() const { return *volfrac; }

    static void variableSetUp();

    static void variableCleanUp();
    
protected:
    static void read_params();

    // call class Amr function to get AmrLevel at lev
    NC &getLevel(int lev) { return dynamic_cast<NC &>(parent->getLevel(lev)); }

    // interpolate from fine to coarse
    // and update temprature
    void avgDown();

    // build geometry information
    void buildMetrics();

    amrex::Real estTimeStep();
    
    // Compute initial time step.
    amrex::Real initialTimeStep();

    // core function
    void compute_dSdt(const amrex::MultiFab &S, amrex::MultiFab &dSdt, amrex::Real dt,
                      amrex::EBFluxRegister *fine, amrex::EBFluxRegister *current);

    // print total momentum and energy
    // only used when verbose>=2
    void printTotal() const;

    // geometric properties
    const amrex::MultiFab* volfrac;
    const amrex::MultiCutFab* bndrycent;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> facecent;

    // cell state descriptor
    amrex::iMultiFab level_mask;
    static constexpr int level_mask_interior   = 0; // valid cells
    static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2; // ghost cells not covered
    static constexpr int level_mask_physbnd    = 3; // outside domain

    amrex::EBFluxRegister flux_reg;

    static constexpr int NUM_GROW = 5;

    // conservative states
    enum StateVariable
    {
        Density = 0,
        Xmom,
        Ymom,
        Zmom,
        Eden,
        NUM_STATE
    };

    enum StateType
    {
        State_Type = 0,
        NUM_STATE_TYPE
    };

    static amrex::BCRec phys_bc;

    // Parameters
    static int verbose;
    static amrex::IntVect hydro_tile_size;
    static amrex::Real cfl;
    static int do_reflux;
    static int do_gravity;
    static int refine_max_dengrad_lev;
    static amrex::Real refine_dengrad;
    static amrex::Vector<amrex::RealBox> refine_boxes;
    static std::string time_integration;
    static int refine_cutcells;
};

#endif
