import numpy as np


class Domain:
    def __init__(self, cluster_id: int, dom_id: int, segments):
        self.cluster_id = cluster_id
        self.domain_id = dom_id
        self.segments: np.array = segments
        self.num_segments = 0
        self.num_residues = 0
        self.count_segments_and_residues()
        self.fit_result = None

    def __str__(self):
        return f"Domain ID : {self.domain_id} \n" \
               f"Cluster ID : {self.cluster_id} \n" \
               f"Number of Segments : {self.num_segments} \n" \
               f"Segments List : {self.segments} \n" \
               f"Number of Residues : {self.num_residues}\n"

    def __repr__(self):
        return f"Domain ID : {self.domain_id} \n" \
               f"Cluster ID : {self.cluster_id} \n" \
               f"Number of Segments : {self.num_segments} \n" \
               f"Segments List : {self.segments} \n" \
               f"Number of Residues : {self.num_residues}\n"

    def add_segment(self, segment, use_end_index=False):
        """
        Adds segment into the list of segments. Rather than simply adding an array, finds the index where the
        given segment connects with the existing segments of the domain.
        Example:
        Domain segments = [[6, 53] [97, 142]]
        Given segment = [54, 62]
        --------------------After adding---------------------
        Domain segments = [[6, 62] [97, 142]]
        :param segment: The segment to be added
        :param use_end_index:   Determines if end index of the given segment is used to connect to the segment. This is
                                for when the given segment is at the front of the chain
        Example:
        Domain segments = [[7, 34] [54, 67]]
        Given segment = [0, 6]
        --------------------After adding---------------------
        Domain segments = [[0, 34] [54, 67]]
        :return:
        """
        if use_end_index:
            # Get the index where the segment joins. There should only be one index found.
            index = np.where(self.segments == (segment[1] + 1))
            self.segments[index[0][0]][index[1][0]] = segment[0]
            pass
        else:
            # Get the index where the segment joins. There should only be one index found.
            index = np.where(self.segments == (segment[0] - 1))
            self.segments[index[0][0]][index[1][0]] = segment[1]
        self.segments = np.sort(self.segments, axis=0)
        self.num_residues += segment[1] + 1 - segment[0]

    def count_segments_and_residues(self):
        self.num_segments = self.segments.shape[0]
        self.num_residues = sum(self.segments[:, 1] + 1 - self.segments[:, 0])

