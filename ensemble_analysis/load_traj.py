import MDAnalysis as mda


def load_traj_from_path(topology_path: str, trajectory_path: (str, list) = None):
    """
    Builds an MDA universe with the supplied topology and trajectory paths
    To build n universe with just a topology which has atom identifiers, don't pass trajectory_path
    :param topology_path: path to the topology file for building the MDA universe
    :param trajectory_path: optional - path to the trajectory file(s) for building the MDA universe
                accepts either a string (path to traj) or a list of strings (paths to traj)
    :return: MDA Universe object
    """

    # sanitize inputs and make sure users are supplying .pdb, .xyz, etc. paths#
    if trajectory_path is None:
        mda_universe = mda.Universe(topology_path)
    else:
        mda_universe = mda.Universe(topology_path, trajectory_path)
    return mda_universe
