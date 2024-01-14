import gemmi
import numpy as np
from sklearn.neighbors import KDTree

from Domain import Domain


def build_domains(slide_window_residues: list, segments_dict: dict, min_domain: int):
    """
    Creates the domains in each cluster.
    :return: List of Domain objects
    """
    print("============================================")
    print("Building domains")
    domains = []
    break_cluster = False
    # For each cluster and their segments
    for cluster, cluster_segments in segments_dict.items():
        # print(f"Cluster {cluster} = {cluster_segments}")
        # Create a binary connection matrix
        bin_mat = create_bin_conn_mat(slide_window_residues, cluster_segments)
        # print(f"Binary Matrix: ")
        # print(f"{bin_mat}")
        # List of Numpy arrays of the binary connection matrix after reduction
        reduced_mat = row_reduction(bin_mat)
        # print(f"Reduced Matrix: {reduced_mat}")
        domains = create_domains(cluster_segments=cluster_segments, domains=domains,
                                 reduced_mat=reduced_mat, cluster_id=cluster)

        # print("--------------------------------------------")

    print_domains(domains)

    domains_to_remove = []
    for domain in domains:
        if domain.num_residues < min_domain:
            domains = join_domains(domain,  domains)
            domains_to_remove.append(domain.domain_id)

    if len(domains_to_remove) > 0:
        new_domains = remove_domains(domains, domains_to_remove)
        print_domains(new_domains, True)
        return new_domains, break_cluster

    print_domains(domains, True)

    return domains, break_cluster


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
    # Go through each segment.
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
            # Get coordinates of all atoms from every residue in the to-be-check segment
            checked_segment_atoms_coordinates = get_segment_atoms_coordinates(residues[checked_segment_start_index:checked_segment_end_index])
            # Checks the distances between the segments
            for r in range(curr_segment_start_index, curr_segment_end_index):
                residue_atoms_coordinates, has_side_chain = get_residue_atoms_coordinates(residues[r])
                is_connected = check_connection_between_atoms(residue_atoms_coordinates, checked_segment_atoms_coordinates,
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


def create_domains(cluster_segments: list, domains: list, reduced_mat: list, cluster_id: int):
    """
    Create domains based on the segments
    :param cluster_segments:
    :param domains:
    :param reduced_mat:
    :param cluster_id:
    :return:
    """
    new_domains = domains
    # For each row in the reduction matrix
    for rows in range(len(reduced_mat)):
        # Get the segment for the cluster
        segments = np.array([cluster_segments[r] for r in reduced_mat[rows]])
        domain = Domain(cluster_id, len(domains), segments)
        new_domains.append(domain)
    return new_domains


def join_domains(tiny_domain: Domain, domains: list):
    """
    Takes the tiny domain's (smaller than minimum domain size) segments and absorbs it into another domain.
    :param tiny_domain: Domain that is smaller than minimum domain size
    :param domains: List of Domain objects
    :return new_domains: Updated List of Domain objects
    """
    tiny_domain_segments = tiny_domain.segments

    # Goes through each segment of the tiny domain
    for ts in range(tiny_domain_segments.shape[0]):
        # First obtain the indices of the previous and next residues connected to the segments
        prev_connecting_index = tiny_domain_segments[ts][0] - 1
        next_connecting_index = tiny_domain_segments[ts][1] + 1
        if prev_connecting_index == -1:
            for curr_d in domains:
                if next_connecting_index in domains[curr_d].segments:
                    domains[curr_d].add_segment(tiny_domain_segments[ts], use_end_index=True)
                    print("Done")
                    break
            continue
        # Going through the other domains
        for d in range(len(domains)):
            # Ignore any domains that are from the same cluster as the tiny domain. Ignore the identical domain as well.
            if domains[d].cluster_id == tiny_domain.cluster_id or domains[d].domain_id == tiny_domain.domain_id:
                continue

            if prev_connecting_index in domains[d].segments:
                domains[d].add_segment(tiny_domain_segments[ts])
                break

    return domains


def remove_domains(domains: list, domains_to_remove: list):
    """
    Removes domains that smaller than minimum domain size
    :param domains: List of Domains
    :param domains_to_remove: IDs of the to-be-removed domains in the list
    :return: new_domains
    """
    new_domains = []
    for domain in domains:
        if domain.domain_id in domains_to_remove:
            continue
        new_domains.append(domain)

    return new_domains


def get_segment_atoms_coordinates(residues: list):
    """
    Gets the residues' atom coordinates
    :param residues:
    :return:
    """
    atoms_pos = np.array([np.asarray(a.pos.tolist()) for res in residues for a in res])
    return atoms_pos


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


def check_connection_between_atoms(residue_atoms: np.array, segment_atoms: np.array, side_chain=False):
    dist_criteria = 4 if side_chain else 10
    # euclidean, l2, minkowski, p, manhattan, cityblock, l1, chebyshev, infinity
    kdt = KDTree(segment_atoms, metric="euclidean")
    count = kdt.query_radius(residue_atoms, r=dist_criteria, count_only=True)
    # print(f"Count: {count} {True if any(c > 0 for c in count) else False}")
    return True if any(c > 0 for c in count) else False


def print_domains(domains: list, after_removal=False):
    print("After: ") if after_removal else print("Before: ")
    for d in domains:
        print(d)

