from glob import glob

from clustering_manager import ClusterBuilder
from ensemble_analysis import load_traj
from ensemble_analysis import calc_rmsd

centroids_base_path = 'example_structures/1zvh/centroids/'
structures_base_path = 'example_structures/1zvh/'


def run_clustering_example():

    results_dict = {}

    centroids_path = glob(centroids_base_path + '*.pdb')
    structures_path = glob(structures_base_path + '*.pdb')

    centroids_traj = load_traj.load_traj_from_path(centroids_path[0], centroids_path)

    cluster_example = ClusterBuilder.ClusterBuilder(centroids_path, 1.3)

    for structure_path in structures_path:
        structure = load_traj.load_traj_from_path(structure_path)

        result = calc_rmsd.calc_rmsd(centroids_traj,
                                     structure,
                                     analysis_ranges=['backbone '
                                                      'and ((resnum 26-32) '
                                                      'or (resnum 52-56) '
                                                      'or (resnum 95-102)) '])
        dataset = (structure_path, result[:, 3])

        new_centroid, index = cluster_example.update_centroids(dataset)

        if new_centroid:
            print(f'Threshold for new centroid met, saving {structure_path} as a new centroid at index {index}')
            centroid_to_save = structure.select_atoms("protein")
            centroid_to_save.write(f'example_structures/1zvh/centroids/voronoi_centroid_1zvh_{index:03d}.pdb')
            centroids_path = glob(centroids_base_path + '*.pdb')
            centroids_traj = load_traj.load_traj_from_path(centroids_path[0], centroids_path)
        else:
            print(f'Threshold for new centroid not met,'
                  f'adding {structure_path} to the population of the centroid at index {index}')

        results_dict[structure_path] = index

    print(f'Number of centroids: {cluster_example.n_centroids}')
    print(f'Population of centroids: {cluster_example.centroid_pops}')

    return results_dict


if __name__ == '__main__':
    results_dict = run_clustering_example()
    print(results_dict)