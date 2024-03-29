
#include <NC.H>
#include <NC_F.H>

using namespace amrex;

static Box the_same_box (const Box& b) { return b; }
//static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

using BndryFunc = StateDescriptor::BndryFunc;

// TODO: Add fix temperature wall boudary
//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall(adiabatic)
//
static int scalar_bc[] =
{
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
};

static int norm_vel_bc[] =
{
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd
};

static int tang_vel_bc[] =
{
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd
};

static void set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static void set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);

    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);

    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);

    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);

    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

void
NC::variableSetUp ()
{
    read_params();

    // fortran get parmparse
    init_fort();

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    // use cell_cons_interp
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,NUM_GROW,NUM_STATE,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);

    Vector<BCRec>       bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);
    BCRec bc;
    int cnt = 0;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "E";

    desc_lst.setComponent(State_Type,
                          Density,
                          name,
                          bcs,
                          BndryFunc(nc_denfill, nc_hypfill));

    StateDescriptor::setBndryFuncThreadSafety(true);

    // DEFINE DERIVED QUANTITIES

    derive_lst.add("T",IndexType::TheCellType(),1,
                   nc_dertemp,the_same_box);
    derive_lst.addComponent("T",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("T",desc_lst,State_Type,Xmom,1);
    derive_lst.addComponent("T",desc_lst,State_Type,Ymom,1);
    derive_lst.addComponent("T",desc_lst,State_Type,Zmom,1);
    derive_lst.addComponent("T",desc_lst,State_Type,Eden,1);

    derive_lst.add("p",IndexType::TheCellType(),1,
                   nc_derpres,the_same_box);
    derive_lst.addComponent("p",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("p",desc_lst,State_Type,Xmom,1);
    derive_lst.addComponent("p",desc_lst,State_Type,Ymom,1);
    derive_lst.addComponent("p",desc_lst,State_Type,Zmom,1);
    derive_lst.addComponent("p",desc_lst,State_Type,Eden,1);

    // Velocities
    // get velocity by momentum/density
    derive_lst.add("ux",IndexType::TheCellType(),1,
                   nc_dervel,the_same_box);
    derive_lst.addComponent("ux",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("ux",desc_lst,State_Type,Xmom,1);

    derive_lst.add("uy",IndexType::TheCellType(),1,
                   nc_dervel,the_same_box);
    derive_lst.addComponent("uy",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("uy",desc_lst,State_Type,Ymom,1);

    derive_lst.add("uz",IndexType::TheCellType(),1,
                   nc_dervel,the_same_box);
    derive_lst.addComponent("uz",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("uz",desc_lst,State_Type,Zmom,1);
}

void
NC::variableCleanUp ()
{
    desc_lst.clear();
    derive_lst.clear();
}
