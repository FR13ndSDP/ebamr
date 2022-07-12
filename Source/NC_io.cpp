
#include <NC.H>
#include <AMReX_EBMultiFabUtil.H>

using namespace amrex;

void
NC::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (do_reflux && level > 0) {
        flux_reg.define(grids, papa.boxArray(level-1),
                        dmap, papa.DistributionMap(level-1),
                        geom, papa.Geom(level-1),
                        papa.refRatio(level-1), level, NUM_STATE);
    }

    buildMetrics();
}

void
NC::checkPoint (const std::string& dir, std::ostream& os, VisMF::How how, bool dump_old)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

void
NC::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    BL_PROFILE("NC::writePlotFile()");
    AmrLevel::writePlotFile(dir, os, how);
}
