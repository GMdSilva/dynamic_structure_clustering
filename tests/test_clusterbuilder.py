# test_calc_rmsd.py

import unittest
from glob import glob

from ensemble_analysis import calc_rmsd
from ensemble_analysis import load_traj
from clustering_manager import ClusterBuilder


class TestClusterBuilder(unittest.TestCase):
    def setUp(self) -> None:
        self.initial_centroids = glob('../example_structures/1zvh/*.pdb')
        self.centroid_candidate = '../example_structures/1zvh/voronoi_centroid_1zvh_2.pdb'
        self.centroid_candidates = glob('../example_structures/1zvh/*.pdb')

    def test_cluster_init(self):
        cluster_test = ClusterBuilder.ClusterBuilder(self.initial_centroids[:1], 0.5)
        self.assertEqual(len(cluster_test.centroid_pops), 1)
        self.assertEqual(cluster_test.centroid_pops[0], 1)
        self.assertEqual(cluster_test.n_centroids, 1)

    def test_one_initial_centroid(self):
        traj = load_traj.load_traj_from_path(self.initial_centroids[0])
        ref = load_traj.load_traj_from_path(self.centroid_candidate)
        result = calc_rmsd.calc_rmsd(traj, ref)
        dataset = ('test', result[:, 3])
        cluster_test = ClusterBuilder.ClusterBuilder(self.initial_centroids[:1], 0.01)
        new_centroid, index = cluster_test.update_centroids(dataset)

        self.assertEqual(new_centroid, True)
        self.assertEqual(cluster_test.centroids[-1], 'test')
        self.assertEqual(index, 2)

        self.assertEqual(cluster_test.centroid_pops[0], 1)
        self.assertEqual(cluster_test.centroid_pops[1], 1)

        self.assertEqual(cluster_test.n_centroids, 2)
        self.assertEqual(len(cluster_test.centroid_pops), 2)

    def test_many_initial_centroids(self):
        traj = load_traj.load_traj_from_path(self.initial_centroids[0], self.initial_centroids[0:3])
        ref = load_traj.load_traj_from_path(self.centroid_candidate)

        result = calc_rmsd.calc_rmsd(traj, ref)
        dataset = ('test', result[:, 3])

        cluster_test = ClusterBuilder.ClusterBuilder(self.initial_centroids[0:3], 0.000001)
        new_centroid, index = cluster_test.update_centroids(dataset)

        self.assertEqual(new_centroid, True)
        self.assertEqual(cluster_test.centroids[-1], 'test')
        self.assertEqual(index, len(self.initial_centroids[0:3]) + 1)

        self.assertEqual(cluster_test.centroid_pops[0], 1)
        self.assertEqual(cluster_test.centroid_pops[-1], 1)

        self.assertEqual(cluster_test.n_centroids, len(self.initial_centroids[0:3]) + 1)
        self.assertEqual(len(cluster_test.centroid_pops), len(self.initial_centroids[0:3]) + 1)

    def test_population_increase(self):
        traj = load_traj.load_traj_from_path(self.initial_centroids[0], self.initial_centroids)
        cluster_test = ClusterBuilder.ClusterBuilder(self.initial_centroids, 0.01)

        for candidate in self.centroid_candidates:
            ref = load_traj.load_traj_from_path(candidate)
            result = calc_rmsd.calc_rmsd(traj, ref)
            dataset = ('test', result[:, 3])
            new_centroid, index = cluster_test.update_centroids(dataset)

        self.assertEqual(new_centroid, False)
        self.assertNotEquals(cluster_test.centroids[-1], 'test')
        self.assertEqual(index, len(self.initial_centroids))

        for pop in cluster_test.centroid_pops:
            self.assertEqual(pop, 2)

        self.assertEqual(cluster_test.n_centroids, len(self.initial_centroids))
        self.assertEqual(len(cluster_test.centroid_pops), len(self.initial_centroids))


if __name__ == '__main__':
    unittest.main()