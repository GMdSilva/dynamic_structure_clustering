# test_calc_rmsd.py

import unittest
from glob import glob

from ensemble_analysis import calc_rmsd
from ensemble_analysis import load_traj


class TestCalcRMSD(unittest.TestCase):
    def setUp(self) -> None:
        self.top_file = '../example_structures/1zvh/voronoi_centroid_1zvh_2.pdb'
        self.traj_files_list = glob('../example_structures/1zvh/*.pdb')
        self.ref_file = '../example_structures/1zvh/voronoi_centroid_1zvh_2.pdb'

        self.ref = load_traj.load_traj_from_path(self.ref_file)
        self.traj = load_traj.load_traj_from_path(self.top_file, self.traj_files_list)

    def test_default_align_analysis_range(self):
        result = calc_rmsd.calc_rmsd(self.traj, self.ref)
        self.assertEqual(result.shape[1], 4)

    def test_custom_align_range(self):
        align_range = 'backbone'
        result = calc_rmsd.calc_rmsd(self.traj, self.ref, align_range=align_range)
        self.assertEqual(result.shape[1], 4)

    def test_custom_single_analysis_range(self):
        analysis_range = ['backbone']
        result = calc_rmsd.calc_rmsd(self.traj, self.ref, analysis_ranges=analysis_range)
        self.assertEqual(result.shape[1], 4)

    def test_custom_multiple_analysis_ranges(self):
        analysis_ranges = ['backbone', 'name CB']
        result = calc_rmsd.calc_rmsd(self.traj, self.ref, analysis_ranges=analysis_ranges)
        self.assertEqual(result.shape[1], 5)

    def test_custom_align_range_and_multiple_analysis_ranges(self):
        analysis_ranges = ['backbone', 'name CB', 'name N']
        align_range = 'backbone'
        result = calc_rmsd.calc_rmsd(self.traj, self.ref, align_range = align_range, analysis_ranges=analysis_ranges)
        self.assertEqual(result.shape[1], 6)


if __name__ == '__main__':
    unittest.main()