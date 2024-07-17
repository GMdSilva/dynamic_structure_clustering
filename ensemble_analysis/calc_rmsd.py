import numpy as np
from typing import List

import MDAnalysis as mda
import MDAnalysis.analysis.rms


def calc_rmsd(traj: mda.Universe,
              ref_frame: mda.Universe,
              align_range: str = 'name CA',
              analysis_ranges: List[str] = ['name CA']) -> np.array:
    """
    For each frame in an MDA Universe trajectory, calculates the RMSD of the atoms in analysis_ranges
    vs. the same atoms in an MDA Universe reference frame after superimposing the structures based on align_range

    :param traj: MDA Universe containing the trajectory for which to calculate the RMSD for
    :param ref_frame: MDA Universe containing a single frame for superimposing and RMSD calculations
    :param align_range: the atom selection (MDAnalysis logic) to use for superimposing the trajectory to the reference
    :param analysis_ranges: a list of atom selections (MDAnalysis logic) for measuring the RMSD
    :return: rmsd_data: np.array containing the RMSD values for each frame in traj for each analysis range
            First column: frame #
            Second column: timestep (if provided, otherwise defaults to 1 frame 1 ps)
            Third column: align_range RMSD
            Fourth through nth column: analysis_range RMSD (one range per column)
    """
    R = MDAnalysis.analysis.rms.RMSD(traj,
                                     ref_frame,
                                     # reference must be a single frame!
                                     select=align_range,
                                     # superimpose based on atoms on align_range
                                     groupselections=analysis_ranges)
                                     # calculate RMSD of atoms on analysis_range
                                     # one calculation for each element
    R.run()
    return R.rmsd
