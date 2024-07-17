from typing import List
import numpy as np


class ClusterBuilder:
    def __init__(self,
                 initial_centroids: List[str],
                 new_centroid_threshold: float):
        """
        Initializes the clustering space with the base centroid(s) identifier(s) and a threshold
        :param initial_centroids: id(s) of the starting centroid(s) to initialize the clustering space lattice
        :param new_centroid_threshold: value threshold for selecting new centroids
        """

        if not initial_centroids:
            raise ValueError("Please define at least one initial centroid for the clustering space")

        if new_centroid_threshold == 0:
            raise ValueError("Threshold for selecting new centroids must be not be zero")

        self.centroids = initial_centroids
        self.centroid_pops = [1] * len(initial_centroids)

        self.n_centroids = len(self.centroids)
        self.threshold = new_centroid_threshold

    def update_centroids(self, dataset: np.array) -> [int, bool]:
        """
        Updates the clustering space (adds new centroids or assigns to one of current centroids populations)
        :param dataset: the values of the comparison between the current entity vs. each of the centroids in order
        :return: if we added a new centroid, return True and its index,
             otherwise, return False and the index that the current entity was added to the population
        """
        counter, add_new_centroid = 1, False
        for i in range(0, len(dataset[1])):
            if dataset[1][i] > self.threshold:
                counter += 1
            else:
                self.centroid_pops[i] += 1
                break
        if counter > self.n_centroids:
            self.centroids.append(dataset[0])
            self.centroid_pops.append(1)
            self.n_centroids += 1
            add_new_centroid = True
        return add_new_centroid, counter

