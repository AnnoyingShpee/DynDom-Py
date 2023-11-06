import gemmi
import numpy as np
from sklearn.neighbors import KDTree

from Domain import Domain

backbone_atoms = ["N", "CA", "C", "O"]


def build_domains(slide_window_residues: list, segments_dict: dict, min_domain: int):
    """
    Creates the domains in each cluster.
    :return: List of Domain objects
    """
    domains = []
    break_cluster = False
    # print(segments_dict)
    # For each cluster and their segments
    for cluster, cluster_segments in segments_dict.items():
        # Create a binary connection matrix
        print("============================================")
        print(f"Cluster {cluster} = {cluster_segments}")
        bin_mat = create_bin_conn_mat(slide_window_residues, cluster_segments)
        print(f"Binary Matrix: ")
        # print(f"{bin_mat}")
        # List of Numpy arrays
        reduced_mat = row_reduction(bin_mat)
        print(f"Reduced Matrix: {reduced_mat}")
        # domains, unused_segments = create_domains(cluster_segments=cluster_segments, min_dom_size=min_domain,
        #                                           domains=domains, reduced_mat=reduced_mat, unused=unused_segments)
        domains, min_dom_not_met = create_domains(cluster_segments=cluster_segments, min_dom_size=min_domain,
                                                  domains=domains, reduced_mat=reduced_mat)

        if (not break_cluster) and min_dom_not_met:
            break_cluster = True

        # print(f"Unused = {unused_segments}")
        # print("--------------------------------------------")

    if break_cluster:
        # print(f"Final Domains = {domains}")
        # print(f"Final Unused = {unused_segments}")
        # domains = join_segments(domains, unused_segments)
        domains = join_segments(domains)

    return domains, break_cluster


def get_residue_atoms_coordinates(residue: gemmi.Residue):
    """
    Get coordinates of all atoms of specified residue and checks whether residue contains side chains (Residues
    containing more than N, CA, C, and O).
    :param residue:
    :return:
    """
    side_chain = False
    if len(residue) > 4:
        side_chain = True
    atoms_pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in residue])
    return atoms_pos, side_chain


def get_segment_atoms_coordinates(residues: list):
    """
    Gets the residues' atom coordinates
    :param residues:
    :return:
    """
    atoms_pos = np.array([np.asarray(a.pos.tolist()) for res in residues for a in res])
    return atoms_pos


def check_connection_between_atoms(residue_atoms: np.array, segment_atoms: np.array, side_chain=False):
    dist_criteria = 4 if side_chain else 10
    # euclidean, l2, minkowski, p, manhattan, cityblock, l1, chebyshev, infinity
    kdt = KDTree(segment_atoms, metric="euclidean")
    count = kdt.query_radius(residue_atoms, r=dist_criteria, count_only=True)
    # print(f"Count: {count} {True if any(c > 0 for c in count) else False}")
    return True if any(c > 0 for c in count) else False


def create_bin_conn_mat(residues: list, segments: np.array):
    """
    Creates a binary connection matrix that determines which segments are connected with each other.
    :param residues: List of gemmi.Residue objects from the slide window
    :param segments:    Numpy array of arrays [s, e], where s and e are the indexes where the segment starts and end in the
                        residues list respectively
    :return bin_conn_mat: A matrix of size [N, N] where N is the length of the segments list
    """
    # Binary connection matrix that stores which segments are connected
    bin_conn_mat = np.zeros((segments.shape[0], segments.shape[0]), dtype=int)
    curr_segment = 0
    # Go through each segment
    while curr_segment < segments.shape[0]:
        # Get the start and end indices of the current residues segment
        curr_segment_start_index = segments[curr_segment][0]
        curr_segment_end_index = segments[curr_segment][1] + 1
        # Check whether current residues segment is connected with any other residues segment
        for s in range(segments.shape[0]):
            # No need to check distances between identical segments.
            if curr_segment == s:
                bin_conn_mat[curr_segment][s] = 1
                continue
            # Get the start and end indices of the segment to be checked
            checked_segment_start_index = segments[s][0]
            checked_segment_end_index = segments[s][1] + 1
            # Get coordinates of all atoms from every residue in the segment to be checked
            segment_atoms_coordinates = get_segment_atoms_coordinates(residues[checked_segment_start_index:checked_segment_end_index])
            # Checks the distances between the segments
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


# def create_domains(cluster_segments: list, min_dom_size: int, domains: list, reduced_mat: list, unused: list):
#     new_domains = domains
#     new_unused = unused
#     min_dom_not_met = False
#     for rows in reduced_mat:
#         segments = np.array([cluster_segments[r] for r in rows])
#         if check_segments_size(segments, min_dom_size):
#             domain = Domain(len(domains), segments)
#             new_domains.append(domain)
#         else:
#             new_unused += [s for s in segments]
#     return new_domains, new_unused


def create_domains(cluster_segments: list, min_dom_size: int, domains: list, reduced_mat: list):
    new_domains = domains
    min_dom_not_met = False
    # For each row in the reduction matrix
    for rows in range(len(reduced_mat)):
        # Get the segment for the cluster
        segments = np.array([cluster_segments[r] for r in reduced_mat[rows]])
        domain = Domain(rows, len(domains), segments)
        if domain.num_residues < min_dom_size:
            min_dom_not_met = True
        new_domains.append(domain)
    return new_domains, min_dom_not_met


def check_segments_size(segments: np.array, min_domain_size: int):
    # segments_sum = np.sum([s[1] + 1 - s[0] for s in segments])
    segments_sum = sum(segments[:, 1] + 1 - segments[:, 0])
    return True if segments_sum >= min_domain_size else False
    # for segment in segments:
    #     count += segment[1] + 1 - segment[0]
    #     if count >= min_domain_size:
    #         return True
    # return False


def join_segments(domains: list):
    """
    Joins by checking where the segment connects at
    Takes any unused segments and
    :param domains: List of Domain objects
    :return domains: Updated List of Domain objects
    """
    new_domains = domains
    domains_to_remove = []
    for curr_d in range(len(domains)):
        print(f"{domains[curr_d].segments} - {domains[curr_d].num_residues}")
        if domains[curr_d].segments.shape[0] < 2 and domains[curr_d].num_residues < 20:
            print(f"Use {domains[curr_d].segments}")
            domains_to_remove.append(curr_d)
            segment = domains[curr_d].segments[0]
            prev_index = segment[0] - 1
            next_index = segment[1] + 1
            for checked_d in range(len(new_domains)):
                if checked_d == curr_d:
                    continue
                print(f"Min = {np.min(new_domains[checked_d].segments)}")
                print(f"Max = {np.max(new_domains[checked_d].segments)}")
                if np.max(new_domains[checked_d].segments) == prev_index:
                    new_domains[checked_d].add_segment(segment)
                    continue
                if np.min(new_domains[checked_d].segments) == next_index:
                    new_domains[checked_d].add_segment(segment)
                    continue
    for i in reversed(domains_to_remove):
        new_domains.pop(i)
    for i in range(len(new_domains)):
        new_domains[i].domain_id = i
    return new_domains


# def join_segments(domains: list, unused: np.array):
#     """
#     Joins by checking where the segment connects at
#     Takes any unused segments and
#     :param domains: List of Domain objects
#     :param unused: Numpy array of segments
#     :return:
#     """
#     new_domains: list = domains
#     for u in unused:
#         # print("=============================")
#         # print(u)
#         prev_index = u[0] - 1
#         next_index = u[1] + 1
#         min_values = np.array([np.min(d.segments) for d in domains])
#         max_values = np.array([np.max(d.segments) for d in domains])
#         # print(f"Min = {min_values}")
#         # print(f"Max = {max_values}")
#         if prev_index > np.max(max_values):
#             # print("Prev More")
#             new_domains[np.argmax(max_values)].add_segment(u)
#             continue
#         if next_index < np.min(min_values):
#             # print("Next Less")
#             new_domains[np.argmin(min_values)].add_segment(u)
#             continue
#         prev_index_domain, next_index_domain = None, None
#         prev_hit, next_hit = False, False
#         added = False
#         # Joins the unused segment
#         for d in range(len(domains)):
#             # print(f"{d} : ({prev_index}, {next_index})")
#             segments: np.array = domains[d].segments
#             pi, pj = np.where(segments == prev_index)
#             ni, nj = np.where(segments == next_index)
#             # print(f"P = {pi}, {pj}")
#             # print(f"N = {ni}, {nj}")
#             if len(pi) > 0 and len(pj) > 0:
#                 prev_index_domain = d
#                 prev_hit = True
#             if len(ni) > 0 and len(nj) > 0:
#                 next_index_domain = d
#                 next_hit = True
#             if prev_hit and next_hit:
#                 if prev_index_domain == next_index_domain:
#                     # print(f"Domain {prev_index_domain} add {u}")
#                     new_domains[prev_index_domain].add_segment(u)
#                 else:
#                     new_domains[next_index_domain].add_segment(u)
#                 added = True
#                 break
#         if not added:
#             if prev_hit:
#                 new_domains[prev_index_domain].add_segment(u)
#             elif next_hit:
#                 new_domains[next_index_domain].add_segment(u)
#     # print(f"New Domains = {new_domains}")
#     return new_domains

