import numpy as np


class Domain:
    def __init__(self, cluster_id: int, dom_id: int, segments):
        self.cluster_id = cluster_id
        self.domain_id = dom_id
        self.segments: np.array = np.sort(segments, axis=0)
        self.num_segments = 0
        self.num_residues = 0
        self.count_segments_and_residues()
        self.fit_result = None
        self.rmsd = 0

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

    def add_segment(self, segment, add_at_left_side=False):
        """
        Adds segment into the list of segments. Rather than simply adding an array, finds the index where the
        given segment connects with the existing segments of the domain. The given segment will be added on the right
        side of the segment unless use_end_index is True
        Example:
        Domain segments = [[6, 53] [97, 142]]
        Given segment = [54, 62]
        --------------------After adding---------------------
        Domain segments = [[6, 62] [97, 142]]
        ============================================================================================================
        :param segment: The segment to be added
        :param add_at_left_side:   Determines if end index of the given segment is used to connect to the segment. This is
                                for when the given segment is at the start of the chain
        Example:
        Domain segments = [[7, 34] [54, 67]]
        Given segment = [0, 6]
        --------------------After adding---------------------
        Domain segments = [[0, 34] [54, 67]]
        :return:
        """
        if add_at_left_side:
            # Get the index where the segment joins. There should only be one index found.
            index = np.where(self.segments[:, 0] == (segment[1] + 1))
            # print("add_at_left_side", add_at_left_side)
            # print("segment\n", segment)
            # print("segments\n", self.segments)
            # print("index", index)
            # print("\n")
            self.segments[index[0][0]][0] = segment[0]

        else:
            # Get the index where the segment joins. There should only be one index found.
            index = np.where(self.segments[:, 1] == (segment[0] - 1))
            # If the index to find is: 78
            # If the segments are:
            # [ [34 45]
            #   [78 95] ]
            # Then the index will be returned as
            # ( array([0, 0], dtype=int64),
            #   array([1, 0], dtype=int64) )
            self.segments[index[0][0]][1] = segment[1]
        self.segments = np.sort(self.segments, axis=0)
        self.num_residues += segment[1] + 1 - segment[0]

    def count_segments_and_residues(self):
        self.num_segments = self.segments.shape[0]
        self.num_residues = sum(self.segments[:, 1] + 1 - self.segments[:, 0])

