from pysimm import system, lmps, cassandra
from collections import OrderedDict
import os
import re
import glob
import types


def mc_md(gas_sst, fixed_sst=None, mc_props=None, md_props=None, np=None, prefix=None, **kwargs):
    """pysimm.apps.mc_md

    Performs the iterative hybrid Monte-Carlo/Molecular Dynamics (MC/MD) simulations using :class:`~pysimm.lmps` for
    MD and :class:`~pysimm.cassandra` for MC

    Args:
        gas_sst (list of :class:`~pysimm.system.System`) : list items describe a different molecule to be
            inserted by MC
        fixed_sst (:class:`~pysimm.system.System`) : fixed during the MC steps group of atoms (default: None)
        WARNING: FOLLOWING NOT IN STABLE RELEASE; implemented by Brandon Tapia (bctapia@mit.edu; btapia@vt.edu; btapia1018@gmail.com)
        np : number of threads for LAMMPS simulations (default: None)
        prefix : prefix to call MPI version of LAMMPS (default: None if np=None, mpiexec if np!=None)



    Keyword Args:
        mcmd_niter (int) : number of MC-MD iterations (default: 10)
        sim_folder (str): relative path to the folder with all simulation files (default: 'results')
        mc_props (dictionary) : description of  all MC properties needed for simulations (see
            :class:`~pysimm.cassandra.GCMC` and :class:`~pysimm.cassandra.GCMC.props` for details)
        md_props (dictionary):  description of all Molecular Dynamics settings needed for simulations (see
            :class:`~pysimm.lmps.Simulation` and :class:`~pysimm.lmps.MolecularDynamics` for details)

    Returns:
        :class:`~pysimm.system.System`:
            Final state of the simulated system
    """

    nonrig_group_name = 'nonrigid_b'
    rig_group_name = 'rigid_b'
    n_iter = kwargs.get('mcmd_niter', 10)
    sim_folder = kwargs.get('sim_folder', 'results')
    xyz_fname = os.path.join(sim_folder, '{:}.md_out.xyz')
    lmps_fname = os.path.join(sim_folder, '{:}.before_md.lmps')

    # Define whether the simulations should be continued or start from the scratch
    l = 1
    is_restart = kwargs.get('restart')
    if is_restart:
        for f in glob.glob(lmps_fname.format('*')):
            l = max(l, int(re.match('\A\d+', os.path.split(f)[1]).group()))

        to_purge = glob.glob(os.path.join(sim_folder, '{:}.*'.format(l + 1))) + \
                   glob.glob(os.path.join(sim_folder, '{:}.md*'.format(l)))
        for f in to_purge:
            os.remove(f)
    # Creating fixed polymer system
    fs = None
    if fixed_sst:
        if isinstance(fixed_sst, system.System):
            fs = fixed_sst
            fs.wrap()
        else:
            print('Cannot setup the fixed system for the simulations. Skipping this')

    # Set the one-molecule gas systems
    gases = []
    if gas_sst:
        if isinstance(gas_sst, system.System):
            gases = [gas_sst]
        elif isinstance(gas_sst, types.ListType):
            for g in cassandra.make_iterable(gas_sst):
                if isinstance(g, system.System):
                    gases.append(g)

    if not gases:
        print('There are no gas molecules were specified correctely\nThe gas molecules are needed to start the '
              'MC-MD simulations\nExiting...')
        exit(1)

    css = cassandra.Cassandra(fixed_sst)
    # Set the Monte-Carlo properties:
    mcp = mc_props
    if mcp:
        CHEM_POT = cassandra.make_iterable(mcp.get('Chemical_Potential_Info'))
        if not CHEM_POT:
            print('Missing chemical potential info\nExiting...')
            exit(1)
    else:
        print('Missing the MC Simulation settings\nExiting...')
        exit(1)
    mcp['Start_Type'] = OrderedDict([('species', [1] + [0] * len(CHEM_POT))])

    # Set the Molecular-Dynamics properties:
    sim = None
    mdp = md_props
    if not mdp:
        print('Missing the MD Simulation settings\nExiting...')
        exit(1)

    # De-synchronizing type names of the framework and the gases to avoid consolidation of types that PySIMM system does
    for gi, g in enumerate(gases):
        for pt in g.particle_types:
            pt.name += '_g' + str(gi + 1)

    while l < n_iter + 1:
        # >>> MC (CASSANDRA) step:
        mcp['Run_Name'] = str(l) + '.gcmc'

        css.add_gcmc(species=gases, is_new=True, chem_pot=CHEM_POT,
                     is_rigid=mcp.get('rigid_type') or [False] * len(gases),
                     out_folder=sim_folder, props_file=str(l) + '.gcmc_props.inp', **mcp)

        if is_restart:
            # Set gas particles positions from the .chk file, and update some properties
            css.run_queue[-1].upd_simulation()
            css.system = css.run_queue[-1].tot_sst.copy()
            # Set frame particles position and box size dimension from the .lmps file
            tmp_sst = system.read_lammps(lmps_fname.format(l))
            for p in css.system.particles:
                p.x = tmp_sst.particles[p.tag].x
                p.y = tmp_sst.particles[p.tag].y
                p.z = tmp_sst.particles[p.tag].z
            css.system.dim = tmp_sst.dim
            is_restart = False
        else:
            css.run()
            css.system.write_lammps(lmps_fname.format(l))

        nm_treads = '1'
        if 'OMP_NUM_THREADS' in os.environ.keys():
            nm_treads = os.environ['OMP_NUM_THREADS']
        os.environ['OMP_NUM_THREADS'] = '1'

        # >>> MD (LAMMPS) step:
        sim_sst = css.system.copy()
        sim_sst.write_lammps(os.path.join(sim_folder, str(l) + '.before_md.lmps'))
        sim = lmps.Simulation(sim_sst, print_to_screen=mdp.get('print_to_screen', False),
                              log=os.path.join(sim_folder, str(l) + '.md.log'))

        sim.add(lmps.Init(cutoff=mdp.get('cutoff'),
                          special_bonds=mdp.get('special_bonds'),
                          pair_modify=mdp.get('pair_modify')))

        # custom definitions for the neighbour list updates
        sim.add_custom('neighbor 2.0 bin \nneigh_modify delay 0 every 1 check yes \nrun_style verlet')

        # adding group definitions to separate rigid and non-rigid bodies
        sim.add(lmps.Group('matrix', 'id', css.run_queue[0].group_by_id('matrix')[0]))
        sim.add(lmps.Group(nonrig_group_name, 'id', css.run_queue[0].group_by_id('nonrigid')[0]))
        rigid_mols = css.run_queue[0].group_by_id('rigid')[0]
        if rigid_mols:
            sim.add(lmps.Group(rig_group_name, 'id', rigid_mols))

        # create the description of the molecular dynamics simulation
        if type(mdp.get('timestep')) == list:

            sim.add(lmps.OutputSettings(thermo=mdp.get('thermo'),
                                        dump={'filename': os.path.join(sim_folder, str(l) + '.md.dump'),
                                              'freq': int(mdp.get('dump'))}))

            for it, (t, lng) in enumerate(zip(mdp.get('timestep'), mdp.get('length'))):

                sim.add(lmps.Velocity(style='create'))
                # adding "run 0" line before velocities rescale for correct temperature init of the
                # system with rigid molecules
                if rigid_mols:
                    sim.add_custom('run 0')
                    sim.add(lmps.Velocity(style='scale'))
                # create the NVT fix for rigid molecules that cannot be put in NPT fix
                if rigid_mols:
                    sim.add(lmps.MolecularDynamics(name='rig_fix_{}'.format(it),
                                                   ensemble='rigid/nvt/small molecule',
                                                   timestep=t,
                                                   length=mdp.get('length'),
                                                   group=rig_group_name,
                                                   temperature=mdp.get('temp'),
                                                   pressure=mdp.get('pressure'),
                                                   run=False))
                sim.add_md(lmps.MolecularDynamics(name='main_fix_{}'.format(it),
                                                  group=nonrig_group_name if rigid_mols else 'all',
                                                  ensemble='npt',
                                                  timestep=t,
                                                  temperature=mdp.get('temp'),
                                                  pressure=mdp.get('pressure'),
                                                  run=False,
                                                  extra_keywords={'dilate': 'all'} if rigid_mols else {}))
                
                sim.add_custom('fix        tether_fix_{} matrix spring tether 30.0 0.0 0.0 0.0 0.0'.format(it))

                sim.add_custom('run {:}\n'.format(lng))
                sim.add_custom('unfix main_fix_{:}'.format(it))
                sim.add_custom('unfix rig_fix_{:}'.format(it))
                sim.add_custom('unfix tether_fix_{:}'.format(it))

        else:
            # create the NVT fix for rigid molecules that cannot be put in NPT fix
            if rigid_mols:
                sim.add(lmps.MolecularDynamics(name='rig_fix',
                                               ensemble='rigid/nvt/small molecule',
                                               timestep=mdp.get('timestep'),
                                               length=mdp.get('length'),
                                               group=rig_group_name,
                                               temperature=mdp.get('temp'),
                                               pressure=mdp.get('pressure'),
                                               run=False))
            sim.add_md(lmps.MolecularDynamics(name='main_fix',
                                              group=nonrig_group_name if rigid_mols else 'all',
                                              ensemble='npt',
                                              timestep=mdp.get('timestep'),
                                              temperature=mdp.get('temp'),
                                              pressure=mdp.get('pressure'),
                                              run=False,
                                              extra_keywords={'dilate': 'all'} if rigid_mols else {}))

            # add the "spring tether" fix to the geometrical center of the system to avoid system creep
            sim.add_custom('fix tether_fix matrix spring tether 30.0 0.0 0.0 0.0 0.0')
            sim.add(lmps.OutputSettings(thermo=mdp.get('thermo'),
                                        dump={'filename': os.path.join(sim_folder, str(l) + '.md.dump'),
                                              'freq': int(mdp.get('dump'))}))
            sim.add_custom('run {:}\n'.format(mdp.get('length')))

        # The input for correct simulations is set, starting LAMMPS:
        if np == None:
            print('[INFO]: np not specified, defaulting to serial')
            sim.run(prefix='[]')
        else:
            if prefix == None:
                print('[INFO]: Prefix not specified, defaulting to mpiexec')
                prefix = 'mpiexec'
            sim.run(np=np, prefix=prefix)
        os.environ['OMP_NUM_THREADS'] = nm_treads

        # Updating the size of the fixed system from the MD simulations and saving the coordinates for the next MC
        # css.system.dim = sim.system.dim
        css.system = sim.system.copy()
        css.unwrap_gas()
        css.system.write_xyz(xyz_fname.format(l))

        mcp['Start_Type']['file_name'] = xyz_fname.format(l)
        mcp['Start_Type']['species'] = [1] + css.run_queue[-1].mc_sst.made_ins
        l += 1

    return sim.system if sim else None

