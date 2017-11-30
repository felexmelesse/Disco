from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py as h5
import tracerUtil as tu
import discoUtil as du
import argparse as ag


def buildFile( xyzFile, repFile , outFile ):
    #In future, allow user to specify which parameters to include and not
    #include in new file to be output

    f = open( xyzFile, "r" )
    g = open( repFile, "r" )
    h = open( outFile, "w")
    step = 0

    line = "start"
    while line !="":
        line = f.readline()
        if line == "":
            break
        Ntr = int( line )
        skip = f.readline()
        line = "{0:d} \nAtoms. Timestep: {1:d}\n".format(Ntr, step)
        h.write( line )
        for i in range( Ntr ):
            line  = f.readline()
            step, time, Id, Type, xxx, yyy, zzz, rrr, phi = tu.getTracerTxtData(line)
            if step==0:
                Lp = 0.0
                Ss = 0.0
            else:
                linea = g.readline()
                repstep, repId, Lp, Ss = tu.getTracerRepData( linea )
                if step != repstep:
                    print "Steps don't match!"
                    print "fStep: ", step
                    print "rStep: ", repstep
                if Id != repId:
                    print "Id's don't match!"
                    print "fId: ", Id
                    print "rId: ", repId
                    quit()
            writeline = "{0:g} {1:d} {2:d} {3:g} {4:g} {5:g} {6:g} {7:g} \n".format(time, step, Id, xxx, yyy, zzz, Lp, Ss )
            h.write( writeline )
        print "Wrote step ", step

    f.close()
    g.close()
    h.close()

    return

if __name__ == '__main__':
    parser = ag.ArgumentParser(description="(re)Write Tracer Output File for Visualization.")
    parser.add_argument('inFile', help="Input file (.xyz) to be editted for visualization.")
    parser.add_argument('repFile', help='Tracer Report file (.dat)')
    parser.add_argument('outFile', help="Output File name")

    args = parser.parse_args()
    inf  = args.inFile
    outf = args.outFile
    repf = args.repFile

    buildFile(inf, repf, outf)
