#ifndef _NC_H_
#define _NC_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

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

    static void variableSetUp();

    static void variableCleanUp();
    
protected:
    static void read_params();

    // call class Amr function to get AmrLevel at lev
    NC &getLevel(int lev) { return dynamic_cast<NC &>(parent->getLevel(lev)); }

    // interpolate from fine to coarse
    // and update temprature
    void avgDown();

    amrex::Real estTimeStep();
    
    // Compute initial time step.
    amrex::Real initialTimeStep();

    // core function
    void compute_dSdt(const amrex::MultiFab &S, amrex::MultiFab &dSdt, amrex::Real dt,
                      amrex::FluxRegister *fine, amrex::FluxRegister *current);

    // print total momentum and energy
    // only used when verbose>=2
    void printTotal() const;

    amrex::FluxRegister flux_reg;

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

    static int refine_max_dengrad_lev;
    static amrex::Real refine_dengrad;

    static amrex::Vector<amrex::RealBox> refine_boxes;
};

#endif
