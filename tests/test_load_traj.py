# test_load_traj.py

import unittest
from glob import glob
import MDAnalysis as mda

from ensemble_analysis import load_traj


class TestTrajLoader(unittest.TestCase):
    def setUp(self):
        self.top_file = '../example_structures/1zvh/voronoi_centroid_1zvh_2.pdb'
        self.traj_files_list = glob('../example_structures/1zvh/*.pdb')

    def test_make_universe_from_pdb_topology_only(self):
        universe = load_traj.load_traj_from_path(self.top_file)
        self.assertIsInstance(universe, mda.Universe)
        self.assertEquals(len(universe.trajectory), 1)

    def test_make_universe_from_topology_and_str_traj(self):
        traj_file_str = self.traj_files_list[0]
        universe = load_traj.load_traj_from_path(self.top_file, traj_file_str)
        self.assertIsInstance(universe, mda.Universe)
        self.assertEquals(len(universe.trajectory), 1)

    def test_make_universe_from_topology_and_list_traj(self):
        universe = load_traj.load_traj_from_path(self.top_file, self.traj_files_list)
        self.assertIsInstance(universe, mda.Universe)
        self.assertEquals(len(universe.trajectory), len(self.traj_files_list))

if __name__ == '__main__':
    unittest.main()

