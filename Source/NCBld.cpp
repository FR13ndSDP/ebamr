
#include <AMReX_LevelBld.H>
#include <NC.H>

using namespace amrex;

class NCBld : public LevelBld
{
    virtual void variableSetUp() override;
    virtual void variableCleanUp() override;
    virtual AmrLevel *operator()() override;
    virtual AmrLevel *operator()(Amr &papa,
                                 int lev,
                                 const Geometry &level_geom,
                                 const BoxArray &ba,
                                 const DistributionMapping &dm,
                                 Real time) override;
};

NCBld nc_bld;

LevelBld*
getLevelBld ()
{
    return &nc_bld;
}

void NCBld::variableSetUp ()
{
    NC::variableSetUp();
}

void
NCBld::variableCleanUp ()
{
    NC::variableCleanUp();
}

AmrLevel*
NCBld::operator() ()
{
    return new NC;
}

AmrLevel*
NCBld::operator() (Amr&            papa,
                         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new NC(papa, lev, level_geom, ba, dm, time);
}
