import gemmi
import numpy as np
from sklearn.neighbors import KDTree

from Domain import Domain

backbone_atoms = ["N", "CA", "C", "O"]


def build_domains(residues: list, segments_dict: dict, min_domain: int):
    """
    Creates the domains
    :return: List of Domain objects
    """
    domains = []
    unused = []
    print(segments_dict)
    for cluster, cluster_segments in segments_dict.items():
        # print(f"{cluster} = {cluster_segments}")
        bin_mat = create_bin_conn_mat(residues, cluster_segments)
        # print("============================================")
        # print("Binary Matrix")
        # print(bin_mat)
        reduced_mat = row_reduction(bin_mat)
        # print("--------------------------------------------")
        # print("Reduced Matrix")
        # print(reduced_mat)
        domains = create_domain(reduced_mat, cluster_segments, min_domain, domains)
    return domains


def get_residue_atoms_coordinates(residue: gemmi.Residue):
    side_chain = False
    if len(residue) > 4:
        side_chain = True
    atoms_pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in residue])
    return atoms_pos, side_chain


def get_segment_atoms_coordinates(residues: list):
    atoms_pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for res in residues for a in res])
    return atoms_pos


def check_connection_between_atoms(residue_atoms: np.array, segment_atoms: np.array, side_chain=False):
    dist_criteria = 4 if side_chain else 10
    # euclidean, l2, minkowski, p, manhattan, cityblock, l1, chebyshev, infinity
    kdt = KDTree(segment_atoms, metric="euclidean")
    count = kdt.query_radius(residue_atoms, r=dist_criteria, count_only=True)
    # print(f"Count: {count} {True if any(c > 0 for c in count) else False}")
    return True if any(c > 0 for c in count) else False


def create_bin_conn_mat(residues: list, segments: list):
    """
    Creates a binary connection matrix that shows which segments are connected with each other.
    :param residues: List of gemmi.Residue objects
    :param segments:    List of tuples (s, e), where s and e are the indexes where the segment starts and end in the
                        residues list respectively
    :return bin_conn_mat: A matrix of size [N, N] where N is the length of the segments list
    """
    # Binary connection matrix that stores which segments are connected
    bin_conn_mat = np.zeros((len(segments), len(segments)), dtype=int)
    curr_segment = 0
    # Go through each segment
    while curr_segment < len(segments):
        # Get the start and end index of the current residues segment
        curr_segment_start_index = segments[curr_segment][0]
        curr_segment_end_index = segments[curr_segment][1] + 1
        # Check whether current residues segment is connected with any other residues segment
        for s in range(len(segments)):
            # No need to check distances between identical segments.
            if curr_segment == s:
                bin_conn_mat[curr_segment][s] = 1
                continue
            # Get the indexes
            checked_segment_start_index = segments[s][0]
            checked_segment_end_index = segments[s][1] + 1
            segment_atoms_coordinates = get_segment_atoms_coordinates(residues[checked_segment_start_index:checked_segment_end_index])
            for r in range(curr_segment_start_index, curr_segment_end_index):
                residue_atoms_coordinates, has_side_chain = get_residue_atoms_coordinates(residues[r])
                is_connected = check_connection_between_atoms(residue_atoms_coordinates, segment_atoms_coordinates,
                                                              has_side_chain)
                if is_connected:
                    bin_conn_mat[curr_segment][s] = 1
                    bin_conn_mat[s][curr_segment] = 1
                    break
        curr_segment += 1
    return bin_conn_mat


# Returns a list
def row_reduction(matrix):
    """
    Performs row reduction on the binary connection matrix to obtain the groups of connected segments. The method uses
    the OR-wise logical operation between the rows containing connections, represented by 1.
    :param matrix: 2D Numpy array
    :return: A list of Numpy arrays containing indices of the segments
    """
    rows = matrix.shape[0]
    connections_list: list = []
    # Go through each row sequentially starting from the row with the highest index to the lowest index
    for m in range(rows-1, -1, -1):
        # The current row will be used to perform OR-wise operations on other rows
        or_wise_result = matrix[m]
        # Go through each index sequentially starting from the current row number to the lowest index
        for n in range(m, -1, -1):
            # OR-wise operations must only be done between rows where the location (m, n) in the matrix is 1
            if matrix[m][n] == 1:
                # Obtain the array at row n
                list_to_or_wise = matrix[n]
                # Perform OR wise logical operation between the 2 arrays
                or_wise_result |= list_to_or_wise
                # Returns a tuple (a, b) where a is the list of indices of or_wise_result where the element is 1.
                # b is empty.
                results = np.where(or_wise_result == 1)
                indices_array = results[0]
                ind = [i for i in range(len(connections_list)) if any(np.isin(indices_array, connections_list[i]))]
                if len(ind) > 0:
                    temp = np.append(connections_list[ind[0]], indices_array)
                    connections_list[ind[0]] = np.unique(temp, return_counts=False)
                    pass
                else:
                    connections_list.append(indices_array)

    return connections_list


def create_domain(reduced_mat: list, cluster_segments: list, min_dom_size: int, domains: list):
    new_domains = domains
    for rows in reduced_mat:
        segments = [cluster_segments[r] for r in rows]
        domain = Domain(len(domains), len(rows), segments)
        if domain.check_domain_size(min_dom_size):
            new_domains.append(domain)
    for i in new_domains:
        print(i)
    return new_domains


# Returns a dictionary
# def row_reduction(matrix):
#     """
#     Performs row reduction on the binary connection matrix to group the connected sets
#     :param matrix: 2D Numpy array
#     :return: A dictionary of connected sets
#     """
#     rows = matrix.shape[0]
#     connections_dict: dict = {}
#     # used_indices = []
#     # Go through each row sequentially starting from the row with the highest index to the lowest index
#     for m in range(rows-1, -1, -1):
#         connected_set = set()
#         lowest_index = m
#         or_wise_result = matrix[m]
#         # Go through each column sequentially starting from the current row number to the lowest index
#         for n in range(m, -1, -1):
#             if matrix[m][n] == 1:
#                 list_to_or_wise = matrix[n]
#                 # Perform OR wise logical operation on the 2 arrays
#                 or_wise_result |= list_to_or_wise
#                 # Returns a tuple (a, b) where a is the list of indices of or_wise_result where the element is 1.
#                 # b is empty.
#                 results = np.where(or_wise_result == 1)
#                 indices_array = results[0]
#                 print(f"Indices = {indices_array}")
#                 lowest_index = indices_array[0]
#                 connected_set.update(set(indices_array))
#                 keys = [key for key, val in connections_dict.items() if any(np.isin(indices_array, list(val)))]
#                 if len(keys) > 0:
#                     connections_dict[min(keys)].update(connected_set)
#                     min_val = min(list(connected_set))
#                     connections_dict[min_val] = connections_dict.pop(min(keys))
#                 else:
#                     # Add a new list of connected segments
#                     connections_dict.update({lowest_index: connected_set})
#                 # print(f"New domain {connections_dict}")
#
#     return connections_dict

