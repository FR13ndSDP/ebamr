
#include <NC.H>

using namespace amrex;

void
NC::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (do_reflux && level > 0) {
        flux_reg.define(grids, dmap, papa.refRatio(level-1), level, NUM_STATE);
    }
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

void NC::printState(const MultiFab& mf)
{
    const IntVect ng(2,2,2);
    const IntVect cell(64,64,0);

    print_state(mf, cell, -1, ng);
}
