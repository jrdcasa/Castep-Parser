import os
import shutil
import datetime
import src.utils.parse_arguments
import src.utils.find_files_recursive
import src.utils.clusterize_rdkit
import src.utils.Logger
from src.cp2k_output import CP2KOutput
from collections import defaultdict


# =============================================================================
def write_summary(energy_dict_sorted, cluster_dict, cluster_centroid_lowest, logger=None):

    hartrees_to_kcalmol = 627.509391

    m = "\t\t           **************** SUMMARY ***************\n"
    print(m) if logger is None else logger.info(m)

    # Write table
    m = "\t\t{0:2s} {1:40s} {2:^19s} {3:^8s} {4:^s}\n".format("#", "Name_job", "Energy(hartrees)",
                                                               "Delta_E (kcal/mol)", "Cluster_Number")
    m1 = "\t\t# "+len(m)*"="
    print(m+m1) if logger is None else logger.info(m+m1)

    m = ""
    e_min = -99999999
    for idx, item in enumerate(energy_dict_sorted.items()):
        key = item[0]
        energy = item[1]
        if idx == 0:
            e_min = energy
        try:
            line = "\t\t{0:^3d} {1:^60s} {2:10.8f} {3:6.2f} {4:^6d}\n".format(idx,
                                                                              key,
                                                                              energy,
                                                                              (energy-e_min)*hartrees_to_kcalmol,
                                                                              cluster_dict[key])
        except (Exception,):
            line = "\t\t{0:^3d} {1:<60s} {2:10.8f} {3:6.2f} {4:s}\n".format(idx,
                                                                            key,
                                                                            energy,
                                                                            (energy-e_min)*hartrees_to_kcalmol,
                                                                            "Not clustering")

        m += line
    print(m) if logger is None else logger.info(m)


# =============================================================================
def write_tcl_pdbs():

    lines = 'proc newRep { sel type color rep imol } {\n'
    lines += '\n'
    lines += '    mol selection $sel\n'
    lines += '    mol representation $type \n'
    lines += '    mol addrep $imol\n'
    lines += '    mol showrep $imol $rep on\n'
    lines += '    mol modcolor $rep $imol $color\n'
    lines += '\n'
    lines += '}\n'
    lines += '\n'
    lines += 'display projection orthographic\n'
    lines += 'axes location on\n'
    lines += 'color Display Background white\n'
    lines += 'display depthcue off\n'
    lines += 'color Labels Bonds black\n'
    lines += 'color Axes Labels black\n'
    lines += '\n'
    lines += 'set listFiles [lsort [glob ./01-PDBs_OUT/*.pdb]]\n'
    lines += 'foreach ifile $listFiles {\n'
    lines += '\n'
    lines += '\n'
    lines += '    mol addfile $ifile  type pdb waitfor -1\n'
    lines += '}\n'
    lines += '    set imol1 [molinfo top]\n'
    lines += 'set nframes [molinfo top get numframes]\n'
    lines += 'mol delrep 0 $imol1\n'
    lines += 'set rep 0\n'
    lines += 'newRep "all and not hydrogen" "DynamicBonds 1.900000 0.100000 12.000000" "name" $rep $imol1 \n'
    lines += 'set rep 1\n'
    lines += 'newRep "all" "VDW 0.200000 12.000000" "name" $rep $imol1 \n'
    lines += 'set rep 2\n'
    lines += 'newRep "all" "DynamicBonds 1.400000 0.100000 12.000000" "name" $rep $imol1 \n'
    lines += 'set rep 3\n'
    lines += 'newRep "all" "HBonds 3.2000 80.0 10.000000" "ColorID 0" $rep $imol1 \n'
    lines += '\n'
    lines += 'animate goto start\n'
    lines += 'pbc box\n'
    lines += '\n'

    with open("VMD_pdb.tcl", 'w') as fvmd:
        fvmd.writelines(lines)


# =============================================================================
def write_gnuplot_templates(energydict, isconvergeddict):

    hartrees_to_kcalmol = 627.509391

    lines = 'reset\n'
    lines += 'set style line 1 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 1.0\n'
    lines += 'set style line 2 lt 1 ps 1.0 lc rgb "red"    pt 5 lw 1.0\n'
    lines += 'set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0\n'
    lines += 'set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0\n'
    lines += 'set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0\n'
    lines += 'set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0\n'
    lines += '\n'
    lines += '###############################################################################\n'
    lines += 'set term wxt 1 enhanced dashed size 1200,800 font "Arial,10"\n'
    lines += 'set multiplot layout 4,4\n'
    lines += '\n'
    lines += 'filelist = system("ls ./02-ITER_DAT/*.data")\n'
    lines += 'set title noenhanced\n'
    lines += 'do for [file in filelist] {\n'
    lines += '  set title sprintf("%s", file)\n'
    lines += '  p file u 1:2 w p ls 1 notitle, file u 1:2 w l ls 1 notitle\n'
    lines += '}\n'
    lines += '\n'
    lines += 'unset multiplot\n'
    lines += '\n'

    with open("iterations.gnu", 'w') as fvmd:
        fvmd.writelines(lines)

    # Write Energy.dat =======================
    with open("./02-ITER_DAT/Energy.dat", 'w') as fvmd:

        lines = '#Structure Label Energy(hartree) DE(kcal/mol) 1/0_Not_converged/COnverged\n'

        emin = None
        for ikey, energy in energydict.items():
            if energy is None:
                continue
            if emin is None:
                emin = energy

        sorted_dict = dict(sorted(energydict.items()))
        for ikey, energy in sorted_dict.items():
            if energy is None:
                continue
            label = ikey.split("_")[-1]
            if isconvergeddict[ikey]:
                cc = 0
            else:
                cc = 1
            lines += '{0:25s} {1:6s} {2:10.5f} {3:10.2f} {4:2d}\n'.format(ikey, label, energy,
                                                                          (energy-emin)*hartrees_to_kcalmol, cc)

        fvmd.writelines(lines)

    lines = 'reset\n'
    lines += 'set style line 1 lt 1 ps 1.0 lc rgb "red"    pt 5 lw 1.0\n'
    lines += 'set style line 2 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 1.0\n'
    lines += 'set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0\n'
    lines += 'set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0\n'
    lines += 'set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0\n'
    lines += 'set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0\n'
    lines += '\n'
    lines += '###############################################################################\n'
    lines += 'set term wxt 1 enhanced dashed size 400,400 font "Arial,8"\n'
    lines += 'set multiplot layout 1,1\n'
    lines += 'set encoding utf8\n'
    lines += '\n'
    lines += 'f1="./02-ITER_DAT/Energy.dat"\n'
    lines += '\n'
    lines += 'set ylabel "{/Symbol D}E_{rel} (kcal/mol)"\n'
    lines += 'set format y  "%.1f"\n'
    lines += 'set key left top\n'
    lines += 'set bmargin at screen 0.25\n'
    lines += '\n'
    lines += '\n'
    lines += 'set xtics 1.0\n'
    lines += 'set mxtics 2\n'
    lines += 'set mytics 2\n'
    lines += 'set xtics rotate by 90 right\n'
    lines += 'set xtics noenhanced\n'
    lines += '\n'
    lines += 'p f1 u ($5==1?$4:$4/0):xtic(1) w p ls 1 title "Not Converged", f1 u 4:xtic(1) w l ls 1 notitle,\\'
    lines += '\n  f1 u ($5==0?$4:$4/0):xtic(1) w p ls 2 title "Converged", f1 u 4:xtic(1) w l ls 2 notitle\n'
    lines += '\n'
    lines += 'unset multiplot\n'
    lines += '\n'

    with open("DeltaEnergy.gnu", 'w') as fvmd:
        fvmd.writelines(lines)


# =============================================================================
def main_app():

    versionapp = 1.0
    # Parse arguments
    args = src.utils.parse_arguments.parse_arguments("CP2K")

    # Setup log
    logger = src.utils.Logger.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Print header
    src.utils.parse_arguments.print_header(versionapp, logger_log=logger)

    energy_dict = defaultdict()
    energy_even_not_converged_dict = defaultdict()
    converged_dict = defaultdict()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t Write PDBs (./01-PDBs_OUT) and DAT (02-ITER_DAT) information {}\n".format(now)
    # Create a directory to store the PDBs
    dir_targets = ["./01-PDBs_OUT", "./02-ITER_DAT", "./03-CLUSTER_XYZ"]
    for idir in dir_targets:
        if os.path.exists(idir):
            shutil.rmtree(idir)
        os.mkdir(idir)

    # Get energy values sorted
    nfiles = 0
    if args.filecp2k is not None:
        c = CP2KOutput(args.filecastep, args.dircp2k, unwrap_coordinates=args.unwrap, logger=logger)
        nfiles = 1
        key = os.path.splitext(os.path.split(args.filecastep)[-1])[0]
        energy_dict[key] = c.getoptenergy()
        converged_dict[key] = c._isconverged
    elif args.dircp2k is not None:
        files = src.utils.find_files_recursive.find_files_matching_pattern(args.dircp2k, args.pattern)
        # Create the objects
        for ifile in files:
            idir = os.path.dirname(os.path.dirname(ifile))
            c = CP2KOutput(ifile, idir, unwrap_coordinates=args.unwrap, logger=logger)
            #print(ifile.replace(args.dircastep, ''))
            #key = os.path.splitext(os.path.split(ifile)[-1])[0]
            key = ifile.replace(idir, '')
            if not c._issuccessopt:
                converged_dict[key] = c._issuccessopt
                try:
                    energy_even_not_converged_dict[key] = c._iteration_energy_hartrees[-1]
                except IndexError:
                    energy_even_not_converged_dict[key] = None
                continue

            energy_dict[key] = c.getoptenergy()
            energy_even_not_converged_dict[key] = c.getoptenergy()
            converged_dict[key] = c._issuccessopt
        nfiles = len(files)
    energy_dict_sorted = dict(sorted(energy_dict.items(), key=lambda x: x[1]))

    m += "\t\t Number of castep files: {}\n".format(nfiles)
    print(m) if logger is None else logger.info(m)

    # Clusterize using the rdkit
    if args.rms_threshold is not None:
        cluster_dict, clustdict_min = src.utils.clusterize_rdkit.clusterize_rdkit(energy_dict_sorted,
                                                                                  args.rms_threshold,
                                                                                  folder="./03-CLUSTER_XYZ",
                                                                                  logger=logger,
                                                                                  cluster_type="RMSD")
    else:
        cluster_dict = defaultdict()
        clustdict_min = defaultdict()
    # Write summary
    write_summary(energy_dict_sorted, cluster_dict, clustdict_min, logger=logger)
    write_tcl_pdbs()
    write_gnuplot_templates(energy_even_not_converged_dict, converged_dict)

    # END JOB
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    version = "0.1"
    main_app()
