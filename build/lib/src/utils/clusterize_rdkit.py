import os
import glob
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdMolAlign, TorsionFingerprints, AllChem
from rdkit.ML.Cluster import Butina


# =============================================================================
def clusterize_rdkit(id_to_energy_dict, threshold=1.5,
                     cluster_type="RMSD", folder="./", logger=None):

    # Type of clustering
    clusters_type_list = ["RMSD", "TFD"]
    if cluster_type.upper() not in clusters_type_list:
        m = '\n\t\t ERROR: Type of clustering "{}" is not allowed.\n'
        m += '\t\t ERROR: Allowing values are RMSD or TFD.\n'
        print(m) if logger is None else logger.error(m)
    cluster_type = cluster_type.upper()

    # Load pdb or mol2
    id_to_file_dict = defaultdict()
    extension = "pdb"
    list_files = sorted(glob.glob(os.path.join("./01-PDBs_OUT", "*.{}".format(extension))))
    nfiles = len(list_files)
    if extension == "mol2":
        molrdkit = Chem.MolFromMol2File(list_files[0], removeHs=False)
    elif extension == "pdb":
        molrdkit = Chem.MolFromPDBFile(list_files[0], removeHs=False)
    else:
        m = '\n\t\t ERROR: "{}" files are not still supported in clusterize rdkit.\n'.format(extension)
        print(m) if logger is None else logger.error(m)
        exit()
    molrdkit.RemoveAllConformers()

    # Assing conformers
    nconf = 0
    filename_to_id = defaultdict()
    for iconf in range(0, nfiles):
        if extension == "mol2":
            mtmp = Chem.MolFromMol2File(list_files[iconf], removeHs=False)
        elif extension == "pdb":
            mtmp = Chem.MolFromPDBFile(list_files[iconf], removeHs=False)
        else:
            m = "\n\t\t ERROR: {} files are not still soported\n".format(extension)
            print(m) if logger is None else logger.error(m)
            exit()
        conf = Chem.Conformer(mtmp.GetConformer(0))
        molrdkit.AddConformer(conf, assignId=True)
        id_to_file_dict[nconf] = list_files[iconf]
        name = os.path.splitext(os.path.split(list_files[iconf])[-1])[0]
        filename_to_id[name] = nconf
        nconf += 1

    # Remove Hydrogens ==================================================================
    m = "\t\t Number of initial conformers: {}\n".format(molrdkit.GetNumConformers())
    m += "\t\t Total Number of atoms: {}\n".format(molrdkit.GetNumAtoms())
    molrdkit_no_h = Chem.RemoveHs(molrdkit)
    m += "\t\t Remove Hydrogens to perform the clustering\n"
    m += "\t\t Heavy Number of atoms: {}\n".format(molrdkit_no_h.GetNumAtoms())
    print(m) if logger is None else logger.info(m)

    # Get all the Best RMS between conformers.
    # GetAllConformerBestRMS does not change the coordinates in molrdkit_no_h
    # The rmsd returned is the lower triangle matrix
    # [ (1,0)
    # [ (2,0) (2,1)
    # [ (3,0) (3,1) (3,2) ...
    m_rmsd_noh_getbest = rdMolAlign.GetAllConformerBestRMS(molrdkit_no_h)

    if cluster_type == "RMSD":
        clusts = Butina.ClusterData(m_rmsd_noh_getbest, nconf, threshold, isDistData=True, reordering=True)
    elif cluster_type == "TFD":
        tfdmat = TorsionFingerprints.GetTFDMatrix(molrdkit)
        clusts = Butina.ClusterData(tfdmat, nconf, threshold, isDistData=True, reordering=True)

    m = "\t\t Number of clusters: {} (threshold: {})\n".format(len(clusts), threshold)
    m += "\t\t Write XYZ clusters and find the lowest energy conformer for each cluster\n"
    print(m) if logger is None else logger.info(m)

    m1 = "\t\t Cluster composition\n"
    ll = int(len(m1))
    m = "\t\t" + ll * "-" + "\n"
    print(m + m1 + m) if logger is None else logger.info(m + m1 + m)

    clustdict_min = defaultdict()
    name_to_cluster = defaultdict()
    # For each cluster
    for i in range(len(clusts)):
        m = "\t\t Cluster {}\n".format(i)
        filecluster = os.path.join(folder, "Cluster{0:03d}.xyz".format(i))
        with open(filecluster, 'w') as fxyz:
            # For each conformer in the cluster
            min_energy = 9999999999999
            min_name = ""
            for jid in clusts[i]:
                name = os.path.split(id_to_file_dict[jid])[-1]
                name_tmp = os.path.splitext(name)[0].split("_")[0:-1]
                name = '_'.join(name_tmp)

                energy = id_to_energy_dict[name]  # in eV
                name_to_cluster[name] = i
                rmsd_noh_centroid = rdMolAlign.GetBestRMS(molrdkit_no_h, molrdkit_no_h, prbId=jid, refId=clusts[i][0])
                rmsd_h_centroid = rdMolAlign.GetBestRMS(molrdkit, molrdkit, prbId=jid, refId=clusts[i][0])
                m += "\t\t\t {0:4d} {1:s} Energy = {2:10.6f} (eV) " \
                     "RMS_from_centroid = {3:6.3f} ({4:6.3f}) (A)\n".format(jid, name,
                                                                            energy, rmsd_noh_centroid,
                                                                            rmsd_h_centroid)
                # Find lowest energy conformer in the cluster
                if energy < min_energy:
                    min_energy = energy
                    min_name = name
                # Write a trajectory xyz file for each cluster to visualize
                pos = molrdkit.GetConformer(jid).GetPositions()
                fxyz.writelines("{}\n".format(len(pos)))
                fxyz.writelines("Cluster{0:03d} - {1:s}\n".format(i, name))
                for idx_atom in range(len(pos)):
                    x, y, z = pos[idx_atom, :]
                    line = "{0:<2s} {1:8.3f} {2:8.3f} {3:10.5f}\n". \
                        format(molrdkit.GetAtomWithIdx(idx_atom).GetSymbol(), x, y, z)
                    fxyz.writelines(line)
            clustdict_min[i] = [min_name, min_energy]
        print(m) if logger is None else logger.info(m)

    # print(AllChem.GetConformerRMS(molrdkit_no_h, 22, 7, prealigned=False))

    # Sort the cluster by energy
    clustdict_min_sorted = sorted(clustdict_min.items(), key=lambda item: item[1][1])

    m = "\t\t           **************** CLUSTER MINIMIZATION ***************\n"
    print(m) if logger is None else logger.info(m)

    # Write table
    m = "\t\t{0:1s} {1:^25s} {2:^19s} {3:^18s}\n".format("#Cluster", "Name_job", "Delta_E (kcal/mol)", "RMS(A) (with H)")
    m1 = "\t\t# "+len(m)*"="
    print(m+m1) if logger is None else logger.info(m+m1)
    m = ""

    # Calculate the RMS of each centroid of the cluster (lowest energy in each cluster).
    # The lowest energy conformer is the reference to calculate the RMS

    ev_to_kcalmol = 23.0609
    eref = clustdict_min_sorted[0][1][1]
    idref = filename_to_id[clustdict_min_sorted[0][1][0]+"_opt"]
    for item in clustdict_min_sorted:
        icluster = item[0]
        name = item[1][0]
        eitem = item[1][1]
        idprobe = filename_to_id[name+"_opt"]
        # Change the structure in the molrdkit. Keep the alligned one.
        rms_noh = rdMolAlign.GetBestRMS(molrdkit_no_h, molrdkit_no_h, prbId=idprobe, refId=idref)
        rms_h = rdMolAlign.GetBestRMS(molrdkit, molrdkit, prbId=idprobe, refId=idref)
        # RMS similar to the calculated in VMD after trajectory alligment
        # rms_noh = AllChem.GetConformerRMS(molrdkit_no_h, idprobe, idref)
        # rms_h = AllChem.GetConformerRMS(molrdkit, idprobe, idref)
        # rms_noh = 0.0
        # rms_h = 0.0

        m += "\t\t{0:4d} {1:^30s} {2:^6.2f} {3:^6.3f} ({4:^6.3f})\n".\
            format(icluster, name, (eitem-eref)*ev_to_kcalmol, rms_noh, rms_h)
    print(m) if logger is None else logger.info(m)

    with open("Cluster_min.xyz", 'w') as fxyz:
        for item in clustdict_min_sorted:
            icluster = item[0]
            name = item[1][0]
            idprobe = filename_to_id[name + "_opt"]
            # Write a trajectory xyz file for each cluster to visualize
            pos = molrdkit.GetConformer(idprobe).GetPositions()
            fxyz.writelines("{}\n".format(len(pos)))
            fxyz.writelines("Cluster{0:03d} - {1:s}\n".format(icluster, name))
            for idx_atom in range(len(pos)):
                x, y, z = pos[idx_atom, :]
                line = "{0:<2s} {1:8.3f} {2:8.3f} {3:10.5f}\n". \
                    format(molrdkit.GetAtomWithIdx(idx_atom).GetSymbol(), x, y, z)
                fxyz.writelines(line)

    lines = 'proc newRep { sel type color rep imol size} {\n'
    lines += '\n'
    lines += '    mol selection $sel\n'
    lines += '    mol representation $type $size\n'
    lines += '    mol addrep $imol\n'
    lines += '    mol showrep $imol $rep on\n'
    lines += '    mol modcolor $rep $imol $color\n'
    lines += '\n'
    lines += '}\n'
    lines += '\n'
    lines += 'display projection orthographic\n'
    lines += 'axes location off\n'
    lines += 'color Display Background white\n'
    lines += 'display depthcue off\n'
    lines += 'color Labels Bonds black\n'
    lines += '\n'
    lines += 'mol new  Cluster_min.xyz type xyz waitfor -1\n'
    lines += 'set imol1 [molinfo top]\n'
    lines += 'set nframes [molinfo top get numframes]\n'
    lines += 'mol delrep 0 $imol1\n'
    lines += 'set rep 0\n'
    lines += 'newRep "all" "cpk" "name" $rep $imol1 0.7\n'
    lines += 'set rep 1\n'
    lines += 'newRep "all" "cpk" "name" $rep $imol1 0.7\n'
    lines += 'mol drawframes $imol1 $rep {0:1:1000}\n'
    lines += 'mol showrep $imol1 $rep 0\n'
    lines += 'mol modcolor $rep $imol1 Timestep\n'
    lines += '\n'
    lines += 'animate goto start\n'
    lines += '\n'
    with open("Cluster_min.tcl", 'w') as fvmd:
        fvmd.writelines(lines)

    return name_to_cluster, clustdict_min_sorted
