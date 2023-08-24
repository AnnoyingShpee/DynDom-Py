

class ScrewAxis:
    def __init__(self, coordinates, vectors, angles, unit_vectors):
        self.coordinates = coordinates
        self.vectors = vectors
        self.angles = angles
        self.unit_vectors = vectors

    def determine_screw_axis(self):
        amptr = self.vectors[:][0] * self.unit_vectors[:][0] + self.vectors[:][1] * self.unit_vectors[:][1] + self.vectors[:][2] * self.unit_vectors[:][2]


