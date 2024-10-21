import glob

import numpy as np
import os


# ------------------------------------------------------------------------------
class CP2KOutput(object):

    """
    Class to get information from a .dat file.
    """

    __slots__ = ['_fname', '_logger', '_issuccessopt', '_iteration_energy_hartrees', '_boxlength', '_boxangle',
                 '_natoms', '_unitcell', '_coordinates', '_elements', '_hartrees_to_kcalmol',
                 '_wrap_coordinates']

    # #########################################################################
    def __init__(self, fname, diroutput, writepdb=True, writeiterdat=True, unwrap_coordinates=True, logger=None):

        self._fname = fname
        self._logger = logger
        self._issuccessopt = False
        self._iteration_energy_hartrees = []
        self._boxlength = []
        self._boxangle = []
        self._natoms = 0
        self._unitcell = np.zeros([3, 3], dtype=np.float32)
        self._coordinates = []
        self._elements = []
        self._hartrees_to_kcalmol = 627.509391
        self._wrap_coordinates = not unwrap_coordinates

        # Parse the cp2k file
        self._parse_cp2k_file()

        # Write pdb
        if writepdb:
            key = os.path.splitext(self._fname.replace(diroutput, ''))[0]
            if key[0] == "/":
                key = key[1:]
            key = key.replace("/", ".")
            pat = key+"_opt.pdb"
            path = os.path.join("./01-PDBs_OUT", pat)
            try:
                if self._wrap_coordinates:
                    # Compute transformation matrix for the triclinic box
                    boxvector = self._boxlength + self._boxangle
                    transformation_matrix = self._compute_transformation_matrix(boxvector)
                    positions = np.array(self._coordinates)
                    self._coordinates = self._apply_periodic_bc(positions, transformation_matrix)
                self._write_pdb(filename_pdb=path)
            except (IndexError, ValueError):
                if os.path.isfile(path):
                    os.remove(path)
                self._issuccessopt = False
        # Write iter dat file
        if writeiterdat:
            key = os.path.splitext(self._fname.replace(diroutput, ''))[0]
            if key[0] == "/":
                key = key[1:]
            key = key.replace("/", ".")
            pat = key+".data"
            path = os.path.join("./02-ITER_DAT", pat)
            try:
                self._write_iterinfo(filename_iter=path)
            except IndexError:
                os.remove(path)
                self._issuccessopt = False

    # #########################################################################
    def _parse_cp2k_file(self):

        iteration_label = "Total Energy               ="
        energy_label = "ENERGY"  # This can appear more than once time per iteration
        iteration_lines = []
        optsuccess_label = "GEOMETRY OPTIMIZATION COMPLETED"
        optsuccess_label2 = "run CONVERGED!"
        numberatoms_label = "Atoms:"

        # Open and read file
        #   1. Extract energy for each iteration
        #   2. Get lattice parameters
        #   3. Check if the optimization is successful.
        tmp_energy = []
        with open(self._fname, 'r') as finp:
            lines = finp.readlines()
            for idx, iline in enumerate(lines):
                # Check if the substring exists in the line
                if iteration_label in iline:
                    iteration_lines.append(iline.strip())
                    self._iteration_energy_hartrees.append(float(iline.split()[-1]))
                if energy_label in iline:
                    tmp_energy.append(float(iline.split()[-1]))
                if optsuccess_label2 in iline:
                    self._issuccessopt = True
                if optsuccess_label in iline or optsuccess_label2 in iline:
                    self._issuccessopt = True

                if numberatoms_label in iline:
                    self._natoms = int(iline.split()[-1])
        # The final energy value is not listed under the label "Total Energy               =".
        # Instead, it appears as "ENERGY". Since the "ENERGY" label can appear multiple
        # times per iteration, it is necessary to manually record the last energy value.
        self._iteration_energy_hartrees.append(tmp_energy[-1])

        # In CP2K the coordinates trajectory are stored in a pdb file
        niter = len(self._iteration_energy_hartrees)
        path = os.path.split(self._fname)[0]
        pdb_file = glob.glob(os.path.join(path, "*.pdb"))
        if len(pdb_file) != 1:
            m = "\n\t\t WARNING: More than one PDB file in the directory.\n"
            m += "\t\t\t{}\n".format(path)
            m += "\t\t  Only the first file is used\n"
            m += "\t\t  {}".format(os.path.split(pdb_file[0])[-1])
            print(m) if self._logger is None else self._logger.info(m)

        # Find the index of the last Step label in the pdb
        idx_tmp = []
        with open(pdb_file[0], 'r') as finp:
            lines = finp.readlines()
            flag_line = "Step"
            for idx, iline in enumerate(lines):
                if flag_line in iline:
                    idx_tmp.append(idx)
        idx_start = idx_tmp[-1]
        self._get_simulation_box(lines, idx_start)
        self._get_optimized_structure(lines, idx_start)

    # #########################################################################
    def _get_simulation_box(self, lines, idx_start):

        idx = idx_start + 1
        a = float(lines[idx].split()[1])   # Angstroms
        b = float(lines[idx].split()[2])
        c = float(lines[idx].split()[3])
        alpha = float(lines[idx].split()[4])  # Degrees
        beta = float(lines[idx].split()[5])
        gamma = float(lines[idx].split()[6])
        self._boxlength = [a, b, c]
        self._boxangle = [alpha, beta, gamma]

        self._unitcell[0, 0] = a
        self._unitcell[1, 1] = b
        self._unitcell[2, 2] = c

    # #########################################################################
    def _get_optimized_structure(self, lines, idx_start):

        idx = idx_start + 2
        for iat in range(0, self._natoms):
            iline = lines[idx]
            atom_name = iline.split()[2]
            self._coordinates.append([float(iline.split()[3]), float(iline.split()[4]), float(iline.split()[5])])
            self._elements.append(atom_name)
            idx += 1

    # #########################################################################
    def _write_pdb(self, filename_pdb="test.pdb"):

        """
        Write a pdb file to check the structure.
        Adapted from MDAnalysis software (https://www.mdanalysis.org/)
        """

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        spacegroup = "P -1"
        zvalue = 1
        lines = ""
        with open(filename_pdb, 'w') as fpdb:
            lines += fmt['REMARK'].format('Extract atoms from CP2K output (J.Ramos)')
            lines += fmt['CRYST1'].format(self._boxlength[0], self._boxlength[1], self._boxlength[2],
                                          self._boxangle[0],
                                          self._boxangle[1],
                                          self._boxangle[2],
                                          spacegroup, zvalue)
            lines += 'MODEL 1\n'

            idx = 0
            while idx < self._natoms:
                resname = "MOL"

                chainid = 'A'
                lines += fmt['ATOM'].format(
                    serial=idx + 1,
                    name=self._elements[idx],  # name = self._atom3d_element[idx]
                    altLoc=" ",
                    resName=resname,
                    chainID=chainid,
                    resSeq=1,
                    iCode=" ",
                    pos=[i for i in self._coordinates[idx]],
                    occupancy=0.0,
                    tempFactor=0.0,
                    segID="    ",
                    element=self._elements[idx])
                idx += 1

            lines += 'ENDMDL\n'
            fpdb.write(lines)

    # #########################################################################
    def _write_iterinfo(self, filename_iter="step.dat"):

        with open(filename_iter, 'w') as fiter:
            lines_iter = "#Iteration E_total(hartrees) DE(kcal/mol)\n"
            eref = self._iteration_energy_hartrees[0]
            for idx, ienergy in enumerate(self._iteration_energy_hartrees):
                lines_iter += "{0:4d} {1:20.10f} {2:6.2f}\n".format(idx, ienergy,
                                                                    (ienergy-eref)*self._hartrees_to_kcalmol)

            fiter.writelines(lines_iter)

    # #########################################################################
    def getoptenergy(self):

        return self._iteration_energy_hartrees[-1]

    # #########################################################################
    @staticmethod
    def _compute_transformation_matrix(box_vectors):
        """
        Compute the transformation matrix for an arbitrary simulation box shape.

        Parameters:
        box_vectors (list): List of box vectors (a, b, c) and angles (alpha, beta, gamma) in degrees.

        Returns:
        numpy.ndarray: 3x3 transformation matrix.
        """
        a, b, c, alpha, beta, gamma = box_vectors

        # Convert angles to radians
        alpha_rad = np.radians(alpha)
        beta_rad = np.radians(beta)
        gamma_rad = np.radians(gamma)

        # Compute box transformation matrix
        cos_alpha = np.cos(alpha_rad)
        cos_beta = np.cos(beta_rad)
        cos_gamma = np.cos(gamma_rad)
        sin_gamma = np.sin(gamma_rad)

        matrix = np.array([
            [a, b * cos_gamma, c * cos_beta],
            [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
            [0, 0, c * np.sqrt(
                1 - cos_alpha ** 2 - cos_beta ** 2 - cos_gamma ** 2 + 2 * cos_alpha * cos_beta * cos_gamma) / sin_gamma]
        ])

        return matrix

    # #########################################################################
    @staticmethod
    def _apply_periodic_bc(positions, transformation_matrix):
        """
        Apply periodic boundary conditions for an arbitrary simulation box shape.

        Parameters:
        positions (numpy.ndarray): N x 3 array of particle positions (x, y, z).
        transformation_matrix (numpy.ndarray): 3x3 transformation matrix for box shape.

        Returns:
        numpy.ndarray: Wrapped particle positions after applying PBC.
        """
        # Transform particle coordinates using the inverse of the transformation matrix
        transformed_positions = np.dot(positions, np.linalg.inv(transformation_matrix))

        # Apply periodic boundary conditions
        wrapped_positions = transformed_positions.copy()

        # Apply PBC for each coordinate axis
        for i in range(3):
            wrapped_positions[:, i] %= 1.0

        # Transform back to original coordinate space
        wrapped_positions = np.dot(wrapped_positions, transformation_matrix)

        return wrapped_positions
