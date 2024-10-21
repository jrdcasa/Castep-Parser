import numpy as np
import os


# ------------------------------------------------------------------------------
class CastepOutput(object):

    """
    Class to get information from a .castep file.
    """

    __slots__ = ['_fname', '_logger', '_issuccessopt', '_iteration_energy_ev', '_boxlength', '_boxangle',
                 '_natoms', '_unitcell', '_coordinates', '_elements', '_ev_to_kcalmol',
                 '_wrap_coordinates', '_map_pdb_castep']

    # #########################################################################
    def __init__(self, fname, diroutput, writepdb=True, writeiterdat=True, unwrap_coordinates=True, logger=None):

        self._fname = fname
        self._logger = logger
        self._issuccessopt = False
        self._iteration_energy_ev = []
        self._boxlength = None
        self._boxangle = None
        self._natoms = 0
        self._unitcell = np.zeros([3, 3], dtype=np.float32)
        self._coordinates = []
        self._elements = []
        self._ev_to_kcalmol = 23.0609
        self._wrap_coordinates = not unwrap_coordinates
        self._map_pdb_castep = dict()  # Map the path of the PDB to the key output label

        # Parse the castep file
        self._parse_castep_file()

        # Write pdb
        if writepdb:
            nopath = os.path.split(self._fname)[-1]
            key = os.path.splitext(self._fname.replace(diroutput, ''))[0]
            if key[0] == "/":
                key = key[1:]
            key = key.replace("/", ".")
            pat = key+"_opt.pdb"
            path = os.path.join("./01-PDBs_OUT", pat)
            self._map_pdb_castep[path] = key
            try:
                wrap_coordinates = True
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
            nopath = os.path.split(self._fname)[-1]
            key = os.path.splitext(self._fname.replace(diroutput, ''))[0]
            if key[0] == "/":
                key = key[1:]
            key = key.replace("/", ".")
            pat = key+".dat"
            path = os.path.join("./02-ITER_DAT", pat)
            try:
                self._write_iterinfo(filename_iter=path)
            except IndexError:
                os.remove(path)
                self._issuccessopt = False

    # #########################################################################
    def _parse_castep_file(self):

        iteration_label = "finished iteration"
        iteration_lines = []
        optsuccess_label = "Geometry optimization completed successfully"
        finalconf_label = "Final Configuration:"
        latticeparam_label = "Lattice parameters"
        reallattice_label = "Real Lattice"
        numberatoms_label = "Total number of ions in cell"

        # Open and read file
        #   1. Extract energy for each iteration
        #   2. Get lattice parameters
        #   3. Check if the optimization is successful.
        with open(self._fname, 'r') as finp:
            lines = finp.readlines()
            for idx, iline in enumerate(lines):
                # Check if the substring exists in the line
                if iteration_label in iline:
                    iteration_lines.append(iline.strip())
                    self._iteration_energy_ev.append(float(iline.split()[-2]))
                if optsuccess_label in iline:
                    self._issuccessopt = True
                if reallattice_label in iline:
                    jline = lines[idx+1]
                    kline = lines[idx+2]
                    lline = lines[idx+3]
                    self._unitcell[0, :] = jline.split()[0:3]
                    self._unitcell[1, :] = kline.split()[0:3]
                    self._unitcell[2, :] = lline.split()[0:3]

                if latticeparam_label in iline:
                    a = float(lines[idx+1].split()[2])   # Angstroms
                    b = float(lines[idx+2].split()[2])
                    c = float(lines[idx+3].split()[2])
                    alpha = float(lines[idx+1].split()[5])  # Degrees
                    beta = float(lines[idx+2].split()[5])
                    gamma = float(lines[idx+3].split()[5])
                    self._boxlength = [a, b, c]
                    self._boxangle = [alpha, beta, gamma]
                if numberatoms_label in iline:
                    self._natoms = int(iline.split()[-1])
                if finalconf_label in iline:
                    idx_start = idx
                    self._get_optimized_structure(lines, idx_start)

    # #########################################################################
    def _get_optimized_structure(self, lines, idx_start):

        idx = idx_start + 11
        for iat in range(0, self._natoms):
            iline = lines[idx]
            atom_name = iline.split()[1]
            uvw = [float(iline.split()[3]), float(iline.split()[4]), float(iline.split()[5]) ]
            idx += 1

            omega = np.dot(self._unitcell[0, :], np.cross(self._unitcell[1, :], self._unitcell[2, :]))
            x = self._boxlength[0] * uvw[0] + \
                self._boxlength[1] * np.cos(np.deg2rad(self._boxangle[2])) * uvw[1] + \
                self._boxlength[2] * np.cos(np.deg2rad(self._boxangle[1])) * uvw[2]

            y = self._boxlength[1] * np.sin(np.deg2rad(self._boxangle[2])) * uvw[1] + \
                self._boxlength[2] * ((np.cos(np.deg2rad(self._boxangle[0])) -
                                       np.cos(np.deg2rad(self._boxangle[1])) * np.cos(np.deg2rad(self._boxangle[2]))) /
                                       np.sin(np.deg2rad(self._boxangle[1]))) * uvw[1]

            z = omega / (self._boxlength[0] * self._boxlength[1] * np.sin(np.deg2rad(self._boxangle[2]))) * uvw[2]

            self._coordinates.append([x, y, z])
            self._elements.append(atom_name)

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
        radtodeg = 180 / np.pi
        lines = ""
        with open(filename_pdb, 'w') as fpdb:
            lines += fmt['REMARK'].format('Extract atoms from Castep output (J.Ramos)')
            lines += fmt['CRYST1'].format(self._boxlength[0], self._boxlength[1], self._boxlength[2],
                                          self._boxangle[0],
                                          self._boxangle[1],
                                          self._boxangle[2],
                                          spacegroup, zvalue)
            lines += 'MODEL 1\n'

            idx = 0
            idx_local = 0
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
            lines_iter = "#Iteration E_total(eV) DE(kcal/mol)\n"
            eref = self._iteration_energy_ev[0]
            for idx, ienergy in enumerate(self._iteration_energy_ev):
                lines_iter += "{0:4d} {1:10.5f} {2:6.2f}\n".format(idx, ienergy, (ienergy-eref)*self._ev_to_kcalmol)

            fiter.writelines(lines_iter)

    # #########################################################################
    def getoptenergy(self):

        return self._iteration_energy_ev[-1]

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
