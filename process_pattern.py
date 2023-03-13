#!/usr/bin/env python3
#
# Copyright (C) 2023 NVIDIA Corporation. All rights reserved.
#
# This program processes data patterns from ECTRANS TRGTOL and TRLTOG modules.
# Raw data is obtained from running ECTRANS with data pattern collection patch,
# e.g.,
#   mpirun -np 512 ectrans-benchmark-sp -t 1279 -n 1 -l 137 --vordiv --scders 
#       --uvders --nproma 32 --norms 2>&1 | tee run.log
#
# Author: Yong Qin (yongq@nvidia.com)
# Version: 0.1 (Mar 9, 2023)
#

import getopt
import numpy as np
import re
import sys
from scipy.io import FortranFile


def grep_sarray(module, proc, bufname, filename):
    pattern = f'^TR{module}.*({proc}).*{bufname}_SIZE_\d = (?P<value>[0-9]+)'
    values = []

    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            match = re.match(pattern, line)
            if match:
                value = int(match.group("value"))
                values.append(value)
    return values


def process(module, mode, source, destination):
    """
    Process data pattern file
    Input:
        module:         TRGTOL or TRLTOG module to process
        mode:           send or receive mode to process
        source:         source process
        destination:    destination process
    """
    # run.log
    rfilename = "run.log"
    # pattern file
    pfilename = "tr" + module + "_pattern_" + mode + ".%d" % source

    if ((module == "gtol" and mode == "s") or (module == "ltog" and mode == "r")):
        # 9-column data pattern and column indices
        columns = 9
        # process column
        cproc = 0
        # temp buffer start and end index column
        cbstart = 1
        cbend = 2
        # user buffer name column (integer, resolvable from parray)
        carray = 3
        # user buffer 1st dimension start and end index column
        carray1start = 4
        carray1end = 5
        # user buffer 2nd, 3rd, 4th dimension index column
        carray2 = 6
        carray3 = 7
        carray4 = 8
        # user buffer name (string)
        parray = ["PGP", "PGP2", "PGPUV", "PGP3A", "PGP3B"]
    elif ((module == "gtol" and mode == "r") or (module == "ltog" and mode == "s")):
        # 8-column data pattern and column indices
        columns = 8
        # process column
        cproc = 0
        # temp buffer start and end index column
        cbstart = 1
        cbend = 2
        # temp buffer stride column
        cbstrid = 3
        # user buffer name column (integer, resolvable from parray)
        carray = 4
        # user buffer 1st dimension start and end index column
        carray1start = 5
        carray1end = 6
        # user buffer 2nd dimension index column
        carray2 = 7
        # user buffer name (string)
        parray = ["PGLAT"]

    pfile = FortranFile(pfilename,'r')
    pdata = pfile.read_ints()
    pfile.close()
    pdata = pdata.reshape((len(pdata) // columns, columns))

    print("Summary:")
    print("  pattern file: %s" % pfilename)
    print("  pattern sample line: %s" % pdata[0])
    if mode == "s":
        print("  process %d send to %d destinations: %s" %
            (source, len(np.unique(pdata[:,0])), np.unique(pdata[:,0])))
    else:
        print("  process %d receive from %d destinations: %s" %
            (source, len(np.unique(pdata[:,0])), np.unique(pdata[:,0])))

    print("  user buffers used (shape):")
    for i in np.unique(pdata[:,carray]):
        if columns == 9:
            # size/shape of user buffers read from data pattern file (local) and run.log (global)
            sarray1 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray1end])
            sarray2 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray2])
            sarray3 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray3])
            sarray4 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray4])
            print("    %s [%d, %d, %d, %d] %s" % ((i, parray[i]),
                # local
                sarray1, sarray2, sarray3, sarray4,
                # global
                grep_sarray(module.upper(), source, parray[i].upper(), rfilename)))
        else:
            # size/shape of user buffers read from data pattern file (local) and run.log (global)
            sarray1 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray1end])
            sarray2 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray2])
            print("    %s [%d, %d] %s" % ((i, parray[i]),
                # local
                sarray1, sarray2,
                # global
                grep_sarray(module.upper(), source, parray[i].upper(), rfilename)))

    unique, counts = np.unique(pdata[:,carray1end]-pdata[:,carray1start]+1, return_counts=True)
    print("  user buffer block size and count: %s" % ([(unique[i], counts[i]) for i in range(len(unique))]))

    #print("data source array 2nd-dim: %d %s" % (len(np.unique(pdata[:,6])), np.unique(pdata[:,6])))
    #print("data source array 3rd-dim: %d %s" % (len(np.unique(pdata[:,7])), np.unique(pdata[:,7])))
    #print("data source array 4rd-dim: %d %s" % (len(np.unique(pdata[:,8])), np.unique(pdata[:,8])))

    if destination != 0:
        print()
        print("Detail:")
        if mode == "s":
            print("  data pattern for %d send to %d:" % (source, destination))
        else:
            print("  data pattern for %d receive from %d:" % (source, destination))

        pdata = pdata[np.isin(pdata[:,cproc],destination)]
        print("  user buffers used (shape):")
        for i in np.unique(pdata[:,carray]):
            if columns == 9:
                # size/shape of user buffers used particularly for this destination
                sarray1 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray1end])
                sarray2 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray2])
                sarray3 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray3])
                sarray4 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray4])
                print("    %s [%d, %d, %d, %d]" % ((i, parray[i]), sarray1, sarray2, sarray3, sarray4))
            else:
                # size/shape of user buffers used particularly for this destination
                sarray1 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray1end])
                sarray2 = np.max(pdata[np.isin(pdata[:,carray],i)][:,carray2])
                print("    %s [%d, %d]" % ((i, parray[i]), sarray1, sarray2))

        unique, counts = np.unique(pdata[:,carray1end]-pdata[:,carray1start]+1, return_counts=True)
        print("  user buffer block size and count: %s" % ([(unique[i], counts[i]) for i in range(len(unique))]))
        print(pdata)



def usage(name):
    """
    Usage page
    Input:
        name: application name
    Output:
        None
    """
    print("Usage: " + name + " [options]")
    print("     -d, --detail            detail output")
    print("     -g, --gtol              TRGTOL module")
    print("     -h, --help              help page")
    print("     -l, --ltog              TRLTOG module")
    print("     -p, --procs <P1[,P2]>   source and destination processes")
    print("                                 P1 send to P2 in send mode")
    print("                                 P1 receive from P2 in receive mode")
    print("                                 summary mode when either is zero")
    print("     -r, --receive           source receives from destination")
    print("     -s, --send              source sends to destination")


def main(name, argv):
    # default configurations
    module = "gtol"
    mode = "s"
    source = 401
    destination = 414
    np.set_printoptions(linewidth=np.inf)

    # parse command line options
    try:
        opts, args = getopt.getopt(argv, "dghlp:rs", ["procs="])
    except getopt.GetoptError:
        print("invalid option: ", argv)
        usage(name)
        exit(-1)

    for opt, val in opts:
        if opt in ["-h", "--help"]:
            usage(name)
            exit(0)
        elif opt in ["-d", "--detail"]:
            np.set_printoptions(linewidth=np.inf, threshold=np.inf)
        elif opt in ["-g", "--gtol"]:
            module = "gtol"
        elif opt in ["-l", "--ltog"]:
            module = "ltog"
        elif opt in ["-p", "--procs"]:
            source, destination = list(map(int, val.split(",")))
        elif opt in ["-r", "--receive"]:
            mode = "r"
        elif opt in ["-s", "--send"]:
            mode = "s"

    # sanity check
    if source == 0:
        print("source process cannot be 0")
        exit(-1)
    elif source < 0 or destination < 0 or source > 512 or destination > 512:
        print("source or destination process out of range")
        exit(-1)

    process(module, mode, source, destination)


if __name__ == '__main__': 
    main(sys.argv[0], sys.argv[1:])
