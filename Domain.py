

class Domain:
    def __init__(self, dom_id: int, num_of_segments: int, segments: list):
        self.id = dom_id
        self.num_segments = num_of_segments
        self.segments = segments

    def __str__(self):
        return f"Domain ID : {self.id} \n" \
               f"Number of Segments : {self.num_segments} \n" \
               f"Segments List : {self.segments}"

    def __repr__(self):
        return f"Domain {self.id} with {self.num_segments} segments"

    def check_domain_size(self, min_domain_size):
        count = 0
        for segment in self.segments:
            count += segment[1] + 1 - segment[0]
            if count >= min_domain_size:
                return True
        return False

