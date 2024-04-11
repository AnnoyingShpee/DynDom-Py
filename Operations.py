import math
import numpy as np
import gemmi
import FileMngr
from Clusterer import Clusterer
from Protein import Protein
from scipy.spatial.transform import Rotation
from scipy.stats import multivariate_normal


class Engine:
    def __init__(self, first_pdb: str, second_pdb: str, params=None):
        # Default parameters to be used if not given input
        if params is None:
            self.parameters = {"chain1id": "A",
                               "chain2id": "A",
                               "window": 5,
                               "domain": 20,
                               "ratio": 1.0,
                               "atoms": "backbone"}
        else:
            self.parameters = params
        # Initialise proteins
        self.protein_1: Protein = Protein(first_pdb, self.parameters["chain1id"], self.parameters["atoms"])
        self.protein_2: Protein = Protein(second_pdb, self.parameters["chain2id"], self.parameters["atoms"])
        print(len(self.protein_1.get_polymer()))
        print(len(self.protein_2.get_polymer()))
        self.fitting_protein_1 = None
        self.main_atoms = []
        if self.parameters["atoms"] == "backbone":
            self.main_atoms = ["N", "CA", "C"]
        elif self.parameters["atoms"] == "ca":
            self.main_atoms = ["CA"]
        # A object containing superimposition information between the proteins
        self.chain_superimpose_result = None
        # List of arrays of gemmi.SupResult objects containing superimposition information between each residue
        self.slide_window_superimpose_results = []
        # Array of translation vectors
        self.translation_vecs = None
        # Array of (3 x 3) rotation matrices
        self.rotation_mats = None
        # List of rotation vectors created from rotation_mats
        self.rotation_vecs = None
        # List of unit vectors of the rotation vectors
        self.unit_vectors = None
        # List of angles
        self.angles = None
        # Bending residues indices
        self.bending_residues = set()
        # Sum of all angles
        self.angles_sum = 0.0
        # Scaling of rotation vectors
        self.scaling = 0.0
        # Clusterer object
        self.clusterer = None

    def run(self):
        """
        Runs the entire program. To make it simple, all operations happen based on Protein 1 conforming to Protein 2.
        :return:
        """
        running = True
        while running:
            chain_compatibility = self.check_chain_length()
            # First check whether the chain lengths are similar
            if not chain_compatibility:
                print("Unable to compare sequences. Chains are too different in length")
                return False
            self.get_slide_window_indices(chain_compatibility)
            # Superimposes Protein 2 onto Protein 1, so they are in the same "space". The superimposed protein will be
            # called Protein S.
            self.superimpose_chains()
            # self.print_chains_superimposed_result()
            # Slides a window over Protein 2's residues and Protein S's residues and superimposes them in each window.
            # This gets the superimposition results for each sliding window.
            self.sliding_window_superimpose_residues()
            # Obtain the rotation matrices from the superposition results
            self.get_transformations()
            # Convert rotation matrices to rotation vectors
            self.convert_rot_mats_to_vecs()
            # self.print_rotation_vectors()
            FileMngr.write_rotation_vec_to_pdb(self.protein_1.id, self.protein_1.slide_window_residues,
                                               self.protein_1.slide_window_residues_indices, self.rotation_vecs)
            # self.print_unit_vectors()
            self.clusterer: Clusterer = Clusterer(self.parameters, self.rotation_vecs,
                                                  self.protein_1, self.protein_2, self.main_atoms)
            self.clusterer.cluster()
            print(f"Cluster status = {self.clusterer.clusterer_status}")
            if self.clusterer.clusterer_status == -1:
                self.parameters["window"] = self.parameters["window"] + 2
                continue
            screws = self.determine_domains_screw_axis()
            for screw in screws:
                print("Screw :")
                print(screw)
            self.determine_bending_residues()
            FileMngr.write_w5_info_file(self.protein_1.id, self.protein_2.id,
                                        param=self.parameters,
                                        domains=self.clusterer.domains,
                                        fixed_domain_id=self.clusterer.fixed_domain)
            running = False
        return True

    def check_chain_length(self):
        """
        Checks whether the residues sequence of Protein 1 and Protein 2 are of the same length.
        If they are not but the difference is small, remove some residues.
        Note: Need to find a way to check which residues to remove.
        :return: Boolean on whether the chains are compatible
        """

        # If the protein chains have identical lengths, then nothing else needs to be done.
        if self.protein_1.chain_atoms.shape[0] == self.protein_2.chain_atoms.shape[0]:
            return 1
        elif abs(self.protein_2.chain_atoms.shape[0] - self.protein_1.chain_atoms.shape[0]) <= 5:
            return 2
        else:
            print("Chains are too different. AKA: One of them is too long.")
            return 0

    def get_slide_window_indices(self, chain_compatibility):
        window_size = int(self.parameters["window"])
        window_mid_index = start_index = (window_size - 1) // 2
        final_index = 0
        if chain_compatibility == 1:
            protein_1_residues_length = len(self.protein_1.get_polymer())
            # The index of the middle residue in the last sliding window
            final_index = protein_1_residues_length + 1 - (window_size - window_mid_index)
            self.protein_1.slide_window_residues_indices = (start_index, final_index)
            self.protein_2.slide_window_residues_indices = (start_index, final_index)
        elif chain_compatibility == 2:
            if self.protein_1.chain_atoms.shape[0] < self.protein_2.chain_atoms.shape[0]:
                protein_1_residues_length = len(self.protein_1.get_polymer())
                # The index of the middle residue in the last sliding window
                final_index = protein_1_residues_length + 1 - (window_size - window_mid_index)
            else:
                protein_2_residues_length = len(self.protein_2.get_polymer())
                # The index of the middle residue in the last sliding window
                final_index = protein_2_residues_length + 1 - (window_size - window_mid_index)
            self.protein_1.slide_window_residues_indices = (start_index, final_index)
            self.protein_2.slide_window_residues_indices = (start_index, final_index)
        else:
            print("Something went wrong.")
        slide_window_residues_length = final_index - start_index + 1
        if window_size >= slide_window_residues_length:
            residues_mid_index = (slide_window_residues_length - 1) // 2
            indices = (residues_mid_index, residues_mid_index + 1)
            self.protein_1.slide_window_residues_indices = indices
            self.protein_2.slide_window_residues_indices = indices

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 1 unto Protein 2 using the backbone atoms.
        :return: Protein 2 chain after transformation
        """
        ptype = self.protein_1.residue_span.check_polymer_type()

        fit_1_to_2_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_2.get_polymer(),
                                                                           self.protein_1.get_polymer(),
                                                                           ptype, gemmi.SupSelect.MainChain)
        # print(f"Count = {fit_1_to_2_result.count}")
        # print(f"Center 1 = {fit_1_to_2_result.center1}")
        # print(f"Center 2 = {fit_1_to_2_result.center2}")
        # print(f"RMSD = {fit_1_to_2_result.rmsd}")
        # print(f"Rotation = {fit_1_to_2_result.transform.mat}")
        # print(f"Translation = {fit_1_to_2_result.transform.vec}")
        self.chain_superimpose_result = fit_1_to_2_result
        self.fitting_protein_1: gemmi.ResidueSpan = self.protein_1.get_polymer()
        self.fitting_protein_1.transform_pos_and_adp(fit_1_to_2_result.transform)

    def sliding_window_superimpose_residues(self):
        """
        Slides an array window of a specified size over both protein chains to superimpose the backbone atoms of the
        residues inside the window and getting indices of each residue located in the middle of each window.
        :return:
        """
        window_size = int(self.parameters["window"])
        target_protein_polymer = self.protein_2.get_polymer()  # The protein that Protein S will superimpose onto.
        self.get_utilised_residues(self.protein_2.slide_window_residues_indices)
        residues_length = self.protein_2.slide_window_residues_indices[1] - self.protein_2.slide_window_residues_indices[0] + 1

        # Only used if the window size ends up becoming larger or equal to the length of the chain. Hopefully, this
        # condition isn't met since the length of a protein chain is quite long.
        if residues_length <= window_size:
            target_protein_polymer_atoms_pos = []
            fitting_protein_polymer_atoms_pos = []
            for r in range(len(target_protein_polymer)):
                for a in self.main_atoms:
                    target_protein_polymer_atoms_pos.append(target_protein_polymer[r].sole_atom(a).pos)
                    fitting_protein_polymer_atoms_pos.append(self.fitting_protein_1[r].sole_atom(a).pos)
            self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                   fitting_protein_polymer_atoms_pos))
        else:
            # For each residue in the 2 proteins, get the backbone atoms
            # Instead of going from the first residue to the last residue of the last sliding window, use the already
            # calculated slide_window_residues_indices.
            for r in range(len(self.protein_1.get_polymer()) - window_size + 1):
                target_protein_polymer_atoms_pos = []
                fitting_protein_polymer_atoms_pos = []
                for i in range(r, r+window_size):
                    for a in self.main_atoms:
                        target_protein_polymer_atoms_pos.append(target_protein_polymer[i].sole_atom(a).pos)
                        fitting_protein_polymer_atoms_pos.append(self.fitting_protein_1[i].sole_atom(a).pos)
                # Superimpose and append the result to list
                self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                       fitting_protein_polymer_atoms_pos))

    def get_utilised_residues(self, indices):
        """
        Gets all residues that were used in the middle of each sliding window
        :param indices:
        :return:
        """
        chain_1 = self.protein_1.get_polymer()
        chain_2 = self.protein_2.get_polymer()
        for i in range(indices[0], indices[1]):
            self.protein_1.slide_window_residues.append(chain_1[i])
            self.protein_2.slide_window_residues.append(chain_2[i])

    def get_transformations(self):
        """
        residues_superimpose_results is a list of gemmi.SupResult objects containing information of the superimposition
        between each sliding window of the Proteins. This function is used to extract the numerical data of the objects
        for KMeans clustering. Specifically, the rotation matrix.
        :return:
        """
        # Initialise the numpy arrays
        # Array of rotation matrices
        self.rotation_mats = np.empty(shape=[len(self.slide_window_superimpose_results), 3, 3])
        # Array of translation vectors
        self.translation_vecs = np.empty(shape=[len(self.slide_window_superimpose_results), 3])
        for i in range(len(self.slide_window_superimpose_results)):
            rot_mat: gemmi.Mat33 = self.slide_window_superimpose_results[i].transform.mat
            trans_vec: gemmi.Vec3 = self.slide_window_superimpose_results[i].transform.vec
            self.rotation_mats[i] = np.asarray(rot_mat.tolist())
            self.translation_vecs[i] = np.asarray(trans_vec.tolist())

    def convert_rot_mats_to_vecs(self):
        """
        Convert rotation matrices to rotation vectors, angles in degrees, and unit rotation vectors
        :return:
        """
        self.rotation_vecs = Rotation.from_matrix(self.rotation_mats).as_rotvec(degrees=True)
        self.angles = np.array([np.linalg.norm(i) for i in self.rotation_vecs])
        self.unit_vectors = self.rotation_vecs / np.linalg.norm(self.rotation_vecs)
        self.angles_sum = np.sum(self.angles)
        self.scaling = 20 * (self.angles_sum / (self.protein_1.slide_window_residues_indices[1] - self.protein_1.slide_window_residues_indices[0]))

    def determine_domains_screw_axis(self):
        """
        Calculates the screw axis of the dynamic domains.
        :return:
        """
        # Get the transformations of Protein 2 fitting to Protein 1 relative to the fixed domain
        fixed_domain_r: gemmi.SupResult = self.get_fixed_domain_transformations()
        # Get the slide chain of Protein 1
        original_protein_1_slide_chain: gemmi.ResidueSpan = self.protein_1.get_slide_window_result()
        # Get the slide chain of Protein 2 and apply the transformation calculated from the fixed domain fitting to move
        # Protein 2 to Protein 1's space
        transformed_protein_2_slide_chain: gemmi.ResidueSpan = self.protein_2.get_slide_window_result()
        transformed_protein_2_slide_chain.transform_pos_and_adp(fixed_domain_r.transform)

        self.clusterer.domains[self.clusterer.fixed_domain].rmsd = fixed_domain_r.rmsd

        domain_screw_axes = []

        # Go through each dynamic domain
        for domain in self.clusterer.domains:
            # Skip the fixed domain
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            """
            First thing to do is to fit Protein 1 onto the transformed Protein 2 so that we can find the displacement
            vectors between the initial positions of the dynamic domain atoms of Protein 1 and the final positions of 
            the dynamic domain atoms of Protein 1
            """
            # Prepare an empty Protein chain that will contain the residues of the domain of Protein 1
            original_protein_1_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_1.chain_param)
            # Prepare an empty Protein chain that will contain the residues of the domain of fitted Protein 1 on 2
            transformed_protein_1_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_1.chain_param)
            # Prepare an empty Protein chain that will contain the residues of the domain of fitted Protein 2
            transformed_protein_2_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_2.chain_param)
            # Go through each segment of the current dynamic domain
            for segment in domain.segments:
                # Add the residues
                for i in range(segment[0], segment[1] + 1):
                    original_protein_1_domain_chain.add_residue(original_protein_1_slide_chain[i])
                    transformed_protein_1_domain_chain.add_residue(original_protein_1_slide_chain[i])
                    transformed_protein_2_domain_chain.add_residue(transformed_protein_2_slide_chain[i])

            # print(f"Ori before : {original_protein_1_slide_chain[4].sole_atom('CA').pos}")
            # print(f"Trans before : {transformed_protein_1_domain_chain[4].sole_atom('CA').pos}")

            # r: gemmi.SupResult = gemmi.superpose_positions(original_atoms, transformed_atoms)
            transformed_protein_1_domain_polymer: gemmi.ResidueSpan = transformed_protein_1_domain_chain.get_polymer()
            transformed_protein_2_domain_polymer: gemmi.ResidueSpan = transformed_protein_2_domain_chain.get_polymer()
            ptype = transformed_protein_1_domain_polymer.check_polymer_type()
            # Fit Protein 1 dynamic domain onto Transformed Protein 2 dynamic domain
            r: gemmi.SupResult = gemmi.calculate_superposition(transformed_protein_2_domain_polymer,
                                                               transformed_protein_1_domain_polymer,
                                                               ptype, gemmi.SupSelect.MainChain)
            domain.rmsd = r.rmsd
            # Transform the domain chain
            transformed_protein_1_domain_polymer.transform_pos_and_adp(r.transform)
            # print(f"Ori after : {original_protein_1_slide_chain[4].sole_atom('CA').pos}")
            # print(f"Trans after : {transformed_protein_1_domain_chain[4].sole_atom('CA').pos}")
            # Get the rotation matrix of the domain transformation and convert to rotation vector
            rot_vec = Rotation.from_matrix(np.asarray(r.transform.mat.tolist())).as_rotvec(degrees=True)
            # Get the unit vector of the rotation vector
            unit_rot_vec = rot_vec / math.sqrt(np.sum(rot_vec**2))
            rot_angle = np.linalg.norm(rot_vec)
            # Get the translation vector of the transformation
            trans_vec = np.asarray(r.transform.vec.tolist())

            """
            Assuming the dynamic domains move as a rigid body, the translations for all atoms would be the same. No need
            to use all residues. Just 4 atoms/residues is fine. Either backbone atoms from 4 residues or just 4 atoms.
            """
            original_atom_coords = np.asarray(original_protein_1_domain_chain[0].sole_atom(self.main_atoms[0]).pos.tolist())
            transformed_atom_coords = np.asarray(transformed_protein_1_domain_chain[0].sole_atom(self.main_atoms[0]).pos.tolist())
            disp_vec = transformed_atom_coords - original_atom_coords

            component_value = np.sum(disp_vec * unit_rot_vec)
            parallel_translation = unit_rot_vec * component_value
            # Calculate difference between displacement and parallel translations to get rotational parts
            rotational_part = disp_vec - parallel_translation
            # Calculate the amplitude of rotation
            rotation_amplitude = math.sqrt(np.sum(rotational_part**2))
            # Calculate unit vector in direction of rotational part
            unit_rotational_part = rotational_part/rotation_amplitude

            # Calculate vector in direction from atoms to axis
            cross_prod_axis = np.cross(unit_rot_vec, unit_rotational_part)
            h_tan = 2*math.tan(0.5*rot_angle)
            atoms_to_axis_direction = (rotation_amplitude*cross_prod_axis)/h_tan

            point_on_axis = original_atom_coords + (0.5 * rotational_part) - atoms_to_axis_direction
            # vector_1_1 = point_on_axis - np.array([-0.911, 0.168, -0.376])
            # vector_1_2 = np.array([-0.911, 0.168, -0.376]) - point_on_axis
            # vector_2_1 = point_on_axis - np.array([0.311, 0.309, 0.899])
            # vector_2_2 = np.array([0.311, 0.309, 0.899]) - point_on_axis
            # unit_vector_1_1 = vector_1_1 / math.sqrt(np.sum(vector_1_1**2))
            # unit_vector_1_2 = vector_1_2 / math.sqrt(np.sum(vector_1_2**2))
            # unit_vector_2_1 = vector_2_1 / math.sqrt(np.sum(vector_2_1**2))
            # unit_vector_2_2 = vector_2_2 / math.sqrt(np.sum(vector_2_2**2))
            # print(f"1_1: {unit_vector_1_1}")
            # print(f"1_2: {unit_vector_1_2}")
            # print(f"2_1: {unit_vector_2_1}")
            # print(f"2_2: {unit_vector_2_2}")
            # print(f"Norm point_on_axis = {np.linalg.norm(point_on_axis)}")
            # print(f"Norm screw axis = {np.linalg.norm(rot_vec)}")
            # print(f"Screw axis = {np.sum(screw_axis * np.array([-0.911, 0.168, -0.376]))}\n")
            # print(f"Screw axis = {np.sum(screw_axis * np.array([0.311, 0.309, 0.899]))}\n")
            domain_screw_axes.append((unit_rot_vec, rot_angle, point_on_axis))

            # # How many atoms/residues to use
            # utilise_num = 4
            # # The displacement vectors of atoms at initial and final conformation
            # disp_vecs = np.empty((utilise_num * len(self.main_atoms), 3))
            # # The rotational parts of atoms
            # rot_parts = np.empty((utilise_num * len(self.main_atoms), 3))
            # # The amplitude of the rotational parts
            # amp_rot_parts = np.empty((utilise_num * len(self.main_atoms), 1))
            # # Unit vectors of the rotational parts
            # unit_rot_parts = np.empty((utilise_num * len(self.main_atoms), 3))
            # for i in range(utilise_num):
            #     for a in range(len(self.main_atoms)):
            #         # Calculate the displacement vectors
            #         original_atom_coords = np.asarray(original_protein_1_domain_chain[i].sole_atom(self.main_atoms[a]).pos.tolist())
            #         transformed_atom_coords = np.asarray(transformed_protein_1_domain_chain[i].sole_atom(self.main_atoms[a]).pos.tolist())
            #         disp_vec = transformed_atom_coords - original_atom_coords
            #         disp_vecs[(i * len(self.main_atoms)) + a] = disp_vec
            #         # Obtain the rotational parts by subtracting the translation vectors from the displacement vectors
            #         rot_parts[(i * len(self.main_atoms)) + a] = disp_vec - unit_trans_vec
            #     start_index = (i * len(self.main_atoms))
            #     end_index = (i * len(self.main_atoms)) + len(self.main_atoms)+1
            #     # Calculate the amplitude of the rotational parts.
            #     amp_rot_parts[i] = math.sqrt(np.sum(rot_parts[start_index:end_index]**2))
            #     unit_rot_parts[start_index:end_index] = rot_parts[start_index:end_index] / amp_rot_parts[i]
            #
            # # norm_disp_vecs = disp_vecs / np.linalg.norm(disp_vecs)
            # final_rot_part = np.sum(rot_parts, axis=0) / unit_rot_parts.shape[0]
            # cross_prod_axis = np.cross(final_rot_part, unit_trans_vec)
            # h_tan = 2 * math.tan(rot_angle)
            # avg_amp_rot = np.sum(amp_rot_parts, axis=0) / amp_rot_parts.shape[0]
            # screw_axis = (avg_amp_rot * cross_prod_axis) / h_tan
            # print(f"Rotation angle = {rot_angle}")
            # print(f"Unit Rotation vector = {unit_rot_vec}")
            # # print(f"Unit rotation vector = {unit_rot_parts}")
            # print(f"Cross prod = {cross_prod_axis}")

            # information = (magnitude, rot_angle, screw_axis)

        return domain_screw_axes

    def determine_bending_residues(self):
        fixed_domain = self.clusterer.domains[self.clusterer.fixed_domain]
        fixed_domain_segments = fixed_domain.segments
        fixed_domain_rot_vecs = self.rotation_vecs[fixed_domain_segments[0][0]:fixed_domain_segments[0][1]]

        # Get the rotation vectors of the fixed domain
        for i in range(1, fixed_domain_segments.shape[0]):
            rot_vecs = self.rotation_vecs[fixed_domain_segments[i][0]:fixed_domain_segments[i][1]]
            fixed_domain_rot_vecs = np.append(fixed_domain_rot_vecs, rot_vecs, axis=0)
        # Calculate mean and standard deviation of the fixed domain rotation vectors
        fixed_domain_mean = np.mean(fixed_domain_rot_vecs, axis=0)
        fixed_domain_std = np.std(fixed_domain_rot_vecs)
        fixed_domain_centered_vecs = fixed_domain_rot_vecs - fixed_domain_mean
        fixed_domain_covar = np.cov(fixed_domain_centered_vecs.T)
        fixed_domain_inv_covar = np.linalg.inv(fixed_domain_covar)
        fixed_domain_var = np.diag(fixed_domain_covar)
        # print(f"Fixed Domain Mean = {fixed_domain_mean}")
        # print(f"Fixed Domain STD = {fixed_domain_std}")
        # print(f"Fixed Domain Var = {fixed_domain_var}")
        # print(f"Fixed Domain Covariance = \n{fixed_domain_covar}")
        # print(f"Fixed Domain Inverse Covariance = \n{fixed_domain_inv_covar}")

        # For each dynamic domain,
        for domain in self.clusterer.domains:
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            print(f"Domain {domain.domain_id}")
            dyn_dom_segments = domain.segments
            dyn_dom_rot_vecs = self.rotation_vecs[dyn_dom_segments[0][0]:dyn_dom_segments[0][1]]

            for i in range(1, dyn_dom_segments.shape[0]):
                rot_vecs = self.rotation_vecs[dyn_dom_segments[i][0]:dyn_dom_segments[i][1]]
                dyn_dom_rot_vecs = np.append(dyn_dom_rot_vecs, rot_vecs, axis=0)

            dyn_dom_mean = np.mean(dyn_dom_rot_vecs, axis=0)
            dyn_dom_std = np.std(dyn_dom_rot_vecs)
            dyn_dom_centered_vecs = dyn_dom_rot_vecs - dyn_dom_mean
            dyn_dom_covar = np.cov(dyn_dom_centered_vecs.T)
            dyn_dom_inv_covar = np.linalg.inv(dyn_dom_covar)
            dyn_dom_var = np.diag(dyn_dom_covar)

            # print(f"Dyn Domain Mean = {dyn_dom_mean}")
            # print(f"Dyn Domain STD = {dyn_dom_std}")
            # print(f"Dyn Domain Var = {dyn_dom_var}")
            # print(f"Dyn Domain Covariance = \n{dyn_dom_covar}")
            # print(f"Dyn Domain Inverse Covariance = \n{dyn_dom_inv_covar}")

            # Calculate the indices of the previous and next residues for each segment of the dynamic domain.
            dyn_dom_prev_indices = dyn_dom_segments[:, 0] - 1
            dyn_dom_next_indices = dyn_dom_segments[:, 1] + 1

            # print("Fixed Domain segments: ")
            # print(fixed_domain.segments)
            # print("Dyn Dom segments:")
            # print(dyn_dom_segments)
            # print("Dyn Dom prev indices: ")
            # print(dyn_dom_prev_indices)
            # print("Dyn Dom next indices: ")
            # print(dyn_dom_next_indices)

            # Get the indices of the fixed domain segments that connects the fixed domain to the dynamic domain.
            # 1D Array of booleans where True means next index after fixed domain segment is dyn dom segment.
            fixed_next_is_dyn = np.in1d(fixed_domain.segments[:, 1], dyn_dom_prev_indices)
            fixed_next_is_dyn_ind = np.where(fixed_next_is_dyn)[0]
            # 1D Array of booleans where True means previous index before fixed domain segment is dyn dom segment.
            fixed_prev_is_dyn = np.in1d(fixed_domain.segments[:, 0], dyn_dom_next_indices)
            fixed_prev_is_dyn_ind = np.where(fixed_prev_is_dyn)[0]

            # Get the indices of the dynamic domain segments that connects the dynamic domain to the fixed domain.
            # 1D Array of booleans where True means next index after dyn dom segment is fixed domain segment.
            dyn_next_is_fixed = np.in1d(dyn_dom_next_indices, fixed_domain.segments[:, 0])
            dyn_next_is_fixed_ind = np.where(dyn_next_is_fixed)[0]
            # 1D Array of booleans where True means previous index before dyn dom segment is fixed domain segment.
            dyn_prev_is_fixed = np.in1d(dyn_dom_prev_indices, fixed_domain.segments[:, 1])
            dyn_prev_is_fixed_ind = np.where(dyn_prev_is_fixed)[0]

            print("Fixed domain segment next index is dyn dom segment: ")
            print(fixed_next_is_dyn)
            print(fixed_next_is_dyn_ind)
            print("Fixed domain segment prev index is dyn dom segment: ")
            print(fixed_prev_is_dyn)
            print(fixed_prev_is_dyn_ind)

            print("Dyn dom segment next index is fixed domain segment: ")
            print(dyn_next_is_fixed)
            print(dyn_next_is_fixed_ind)
            print("Dyn dom segment prev index is fixed domain segment: ")
            print(dyn_prev_is_fixed)
            print(dyn_prev_is_fixed_ind)

            p = 4.6

            # Go backwards through the fixed domain residues of the segments
            for segment_ind in fixed_next_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                for i in range(segment[1], segment[1] - 1, -1):
                    centered_vec = self.rotation_vecs[i] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    # print("Backward Fixed Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go forwards through the dyn dom residues of the segments
            for segment_ind in dyn_prev_is_fixed_ind:
                segment = domain.segments[segment_ind]
                for i in range(segment[1], segment[0] + 1):
                    centered_vec = self.rotation_vecs[i] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Forward Dyn Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go forwards through the fixed domain residues of the segments
            for segment_ind in fixed_prev_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                for i in range(segment[0], segment[1] + 1):
                    centered_vec = self.rotation_vecs[i] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    # print("Forward Fixed Q Value =", q_value)
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            # Go backwards through the dyn dom residues of the segments
            for segment_ind in dyn_next_is_fixed_ind:
                segment = domain.segments[segment_ind]
                for i in range(segment[1], segment[0] - 1, -1):
                    centered_vec = self.rotation_vecs[i] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Backward Dyn Q Value =", q_value)
                    # if q_value < p or q_value > 1-p:
                    if q_value > p:
                        self.bending_residues.add(i)
                    else:
                        break

            print(self.bending_residues)
            print(len(self.bending_residues))

    def get_fixed_domain_transformations(self):
        slide_window_1 = self.protein_1.get_slide_window_result()
        slide_window_2 = self.protein_2.get_slide_window_result()
        coords_1 = []
        coords_2 = []
        fixed_domain_id = self.clusterer.fixed_domain
        for s in self.clusterer.domains[fixed_domain_id].segments:
            for i in range(s[0], s[1] + 1):
                for a in self.main_atoms:
                    coords_1.append(slide_window_1[i].sole_atom(a).pos)
                    coords_2.append(slide_window_2[i].sole_atom(a).pos)

        r: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2)
        return r

    def determine_external_component_motion(self):

        return

    def print_chains_superimposed_result(self):
        # A gemmi.SupResult object containing superimposition information between 2 chains
        print(f"RMSD =                  {self.chain_superimpose_result.rmsd}")
        print(f"Count =                 {self.chain_superimpose_result.count}")
        print(f"Center 1 =              {self.chain_superimpose_result.center1}")
        print(f"Center 2 =              {self.chain_superimpose_result.center2}")
        print(f"Translation Vector =    {self.chain_superimpose_result.transform.vec}")
        print(f"Rotation Matrix =       {self.chain_superimpose_result.transform.mat}")

    def print_slide_window_superimpose_results(self, n=None):
        if n is None or n > len(self.slide_window_superimpose_results):
            n = len(self.slide_window_superimpose_results)
        print(f"slide_window_superimpose_result size = {len(self.slide_window_superimpose_results)}")
        for i in range(n):
            item: gemmi.SupResult = self.slide_window_superimpose_results[i]
            print(f"RMSD =                  {item.rmsd}")
            print(f"Count =                 {item.count}")
            print(f"Center 1 =              {item.center1}")
            print(f"Center 2 =              {item.center2}")
            print(f"Translation Vector =    {item.transform.vec}")
            print(f"Rotation Matrix =       {item.transform.mat}")

    def print_slide_window_residue_indices(self):
        print(f"slide_window_residue_indices = {self.protein_1.slide_window_residues_indices}")

    def print_rotation_matrices(self, n=None):
        if n is None or n > self.rotation_mats.shape[0]:
            n = self.rotation_mats.shape[0]
        print(f"rotation_mats shape = {self.rotation_mats.shape}")
        print(f"rotation_mats[0:{n}] = {self.rotation_mats[0:n]}")

    def print_rotation_vectors(self, n=None):
        print(f"rotation_vecs shape = {self.rotation_vecs.shape}")
        if n is None or n > self.rotation_vecs.shape[0]:
            n = self.rotation_vecs.shape[0]
        for i in range(n):
            print(f"{self.protein_1.slide_window_residues[i]} - [{self.rotation_vecs[i][0]}, {self.rotation_vecs[i][1]}, {self.rotation_vecs[i][2]}]")

    def print_unit_vectors(self, n=None):
        print(f"rotation_vecs shape = {self.unit_vectors.shape}")
        if n is None or n > self.unit_vectors.shape[0]:
            n = self.unit_vectors.shape[0]
        for i in range(n):
            print(
                f"{self.protein_1.slide_window_residues[i]} - [{self.unit_vectors[i][0]}, {self.unit_vectors[i][1]}, {self.unit_vectors[i][2]}]")

    # def print_angles(self, n=None):
    #     if n is None or n > self.angles.shape[0]:
    #         n = self.angles.shape[0]
    #     print(f"angles shape = {self.angles.shape}")
    #     print(f"angles[0:{n} = {self.angles[0:n]}")

    #
    # def print_chains_superimposed_result(self):
    #     print(math.sqrt(sum([r.rmsd for r in self.chain_superimpose_result])))
    #     # for r in self.chain_superimpose_result:
    #     #     print(f"RMSD =                  {r.rmsd}")
    #     #     print(f"Count =                 {r.count}")
    #     #     print(f"Center 1 =              {r.center1}")
    #     #     print(f"Center 2 =              {r.center2}")
    #     #     print(f"Translation Vector =    {r.transform.vec}")
    #     #     print(f"Rotation Matrix =       {r.transform.mat}")

