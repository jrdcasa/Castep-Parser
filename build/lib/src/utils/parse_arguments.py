import argparse
import os
import datetime
import sys


# =============================================================================
def parse_arguments(program):

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    list_of_programs = ["CASTEP", "CP2K"]
    if program not in list_of_programs:
        m = '\n\t\tERROR: Program "{}" is not implemented.\n'.format(program.upper())
        for i in list_of_programs:
            m += "\t\t\t{}\n".format(i)
        raise argparse.ArgumentTypeError(m)

    desc = """Parse {} Output files""".format(program)

    parser = argparse.ArgumentParser(description=desc)

    # Create a mutually exclusive group for the two arguments
    group = parser.add_mutually_exclusive_group()

    if program.upper() == "CASTEP":
        group.add_argument("-f", "--file", dest="filecastep",
                           help="A *.castep output file from CASTEP.",
                           action="store", default=None)

        group.add_argument("-d", "--dir", dest="dircastep",
                           help="A directory to search for *.castep output files from CASTEP.",
                           action="store", default="./")

        parser.add_argument("--log", dest="log",
                            help="Name of the file to write logs from this command",
                            action="store", required=False, default="castep_parser.log")
    elif program == "CP2K":
        group.add_argument("-f", "--file", dest="filecp2k",
                           help="A *.dat output file from CP2K.",
                           action="store", default=None)

        group.add_argument("-d", "--dir", dest="dircp2k",
                           help="A directory to search for *.castep output files from CP2K.",
                           action="store", default=None)

        parser.add_argument("-p", "--pattern", dest="pattern",
                            help="A pattern to search for output files from CP2K.",
                            action="store", required=True, default=None)

        parser.add_argument("--log", dest="log",
                            help="Name of the file to write logs from this command",
                            action="store", required=False, default="cp2k_parser.log")

    parser.add_argument("-t", "--threshold", dest="rms_threshold", type=float,
                        help="Threshold rms in angstroms. Default value 1.0 Angstroms.",
                        action="store", default=None)

    parser.add_argument("--unwrap", dest="unwrap",
                        help="Unwrap PDB coordinates.",
                        action="store_true")

    args = parser.parse_args()

    if program.upper() == "CASTEP":
        if args.filecastep is not None and not os.path.isfile(args.filecastep):
            raise argparse.ArgumentTypeError('\nERROR: File "{}" does not exist.\n'.format(args.filecastep))

        if args.filecastep is not None and os.path.splitext(args.filecastep)[-1] != ".castep":
            raise argparse.ArgumentTypeError('\nERROR: Input File "{}" must be a castep file.\n'.format(args.filecastep))

        if args.dircastep is not None and not os.path.isdir(args.dircastep):
            raise argparse.ArgumentTypeError('\nERROR: Directoy "{}" does not exist.\n'.format(args.dircastep))
    elif program.upper() == "CP2K":
        if args.filecp2k is not None and not os.path.isfile(args.filecp2k):
            raise argparse.ArgumentTypeError('\nERROR: File "{}" does not exist.\n'.format(args.filecp2k))

        if args.filecp2k is not None and os.path.splitext(args.filecastep)[-1] != ".dat":
            raise argparse.ArgumentTypeError('\nERROR: Input File "{}" must be a CP2K file.\n'.format(args.filecp2k))

        if args.dircp2k is not None and not os.path.isdir(args.dircp2k):
            raise argparse.ArgumentTypeError('\nERROR: Directoy "{}" does not exist.\n'.format(args.dircp2k))

    return args


# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                               Castep Parser  
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This utility is part of the polyanagro library. Polyanagro is an 
        open-source python library to analyze simulations of polymer systems.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger_log is None else logger_log.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    m += "\t\t\t         or\n"
    m += "\t\t\t{}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)
