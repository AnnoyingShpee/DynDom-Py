import numpy as np


class Domain:
    def __init__(self, dom_id: int, segments):
        self.id = dom_id
        self.segments: np.array = segments
        self.num_segments = 0
        self.num_residues = 0
        self.count_segments_and_residues()

    def __str__(self):
        return f"(Domain ID : {self.id} \n" \
               f"Number of Segments : {self.num_segments} \n" \
               f"Segments List : {self.segments} \n" \
               f"Number of Residues : {self.num_residues})"

    def __repr__(self):
        return f"(Domain ID : {self.id} \n" \
               f"Number of Segments : {self.num_segments} \n" \
               f"Segments List : {self.segments} \n" \
               f"Number of Residues : {self.num_residues})"

    def add_segment(self, segment):
        self.segments = np.append(self.segments, [segment], axis=0)
        self.num_segments += 1
        self.num_residues += segment[1] + 1 - segment[0]

    def count_segments_and_residues(self):
        self.num_segments = self.segments.shape[0]
        self.num_residues = sum(self.segments[:, 1] + 1 - self.segments[:, 0])

