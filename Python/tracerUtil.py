from __future__ import division
import numpy as np
import h5py as h5

def getTracerTxtData( line ):
    arr  = line.split()
    step = int( arr[0] )
    time = float( arr[1] )
    Id   = int(arr[2] )
    Type = int( arr[3] )
    xxx  = float( arr[4] )
    yyy  = float( arr[5] )
    zzz  = float( arr[6] )
    rrr  = float( arr[7] )
    phi  = float( arr[8] )

    return step, time, Id, Type, xxx, yyy, zzz, rrr, phi

def loadTrH5Data( filename ):

    f = h5.File(filename, 'r')

    TrPhaseSp = f['Data']['Tracer_data'][:]
    TrInts     = f['Data']['Tracer_id'][:]

    rpz   = TrPhaseSp[:,:3]
    vvv   = TrPhaseSp[:,3:]
    ids   = TrInts[:,0]
    types = TrInts[:,1]

    f.close()

    return rpz, vvv, ids, types


def find_mini_trs( inFile, stepsize, cutframe ):

    f = open( inFile, "r" )
    r_mini = 1.0

    t_cut = stepsize*cutframe
    step  = 0

    tr_minidisk  = []
    tr_innerdisk = []
    line = "start"
    while line != "":
        line = f.readline()
        if line == "":
            break
        Ntr  = int( line )
        skip = f.readline()
        for i in range( Ntr ):
            line = f.readline()
            step, time, Id, Type, xxx, yyy, zzz, rrr, phi = getTracerData(line)
            Type = 1

            pl1, pl2 = getBinaryPos( time )
            dl2_1 = (rrr - pl1[0])**2 + 2*rrr*pl1[0]*( np.cos(phi - pl1[1]) )**2
            dl2_2 = (rrr - pl2[0])**2 + 2*rrr*pl2[0]*( np.cos(phi - pl2[1]) )**2
            dist1 = np.sqrt( dl2_1 )
            dist2 = np.sqrt( dl2_2 )

            mini   = (dist1<r_mini) or (dist2<r_mini)
            outer  = (step < t_cut)
            if( mini ):
                tr_minidisk.append(Id)
                if( outer ):
                    tr_innerdisk.append(Id)
            print "Finding ", step

    f.close()

    tr_minidisk  = list(set(tr_minidisk))
    print "Total number of tracers in minidisk ", len(tr_minidisk)
    tr_outerdisk = list( set(tr_minidisk) - set(tr_innerdisk) )
    print "Number of tracers from Outer Disk ", len(tr_outerdisk)

    g = open('Id_list.dat', 'w')
    for id in tr_minidisk:
        g.write("%s\n" % id)
    g.close()

    h = open('Id_outlist.dat', 'w')
    for id in tr_outerdisk:
        h.write("%s\n" % id)
    h.close()

    #return tr_minidisk

def tr2mini( inFile, outFile ):

    #tracers = find_mini_trs( inFile )
    step  = 0
    tracers = np.loadtxt("Id_list.dat")
    outers  = np.loadtxt("Id_outlist.dat")
    trnum   = tracers.shape[0] + 1

    f = open( inFile , "r" )
    g = open( outFile, "w" )

    line = "start"
    while line != "":
        line = f.readline()
        if( line=="" ):
            break
        Ntr  = int( line )
        skip = f.readline()
        step = int( skip.split()[2] )
        line = "{0:d} \nAtoms. Timestep: {1:d}\n".format(trnum, step)
        g.write( line )
        for i in range( Ntr ):
            line = f.readline()
            step, time, Id, Type, xxx, yyy, zzz, rrr, phi = getTracerData(line)

            #print time
            Type = 1
            if Id in outers:
                Type = 2

            if Id in tracers:
                line = "{0:d} {1:g} {2:d} {3:d} {4:g} {5:g} {6:g} {7:g} {8:g} \n".format( step, time, Id, Type, xxx, yyy, zzz, rrr, phi)
                g.write( line )
            print "Writing ", step

def reTyping( inFile, outFile ):

    step = 0
    retype = np.loadtxt('Id_outlist.dat')
    outers  = np.loadtxt("Id_outlist.dat")

    f = open( inFile , "r" )
    g = open( outFile, "w" )

    line = "start"
    while line != "":
        line = f.readline()
        if( line=="" ):
            break
        Ntr  = int( line )
        skip = f.readline()
        step = int( skip.split()[2] )
        line = "{0:d} \nAtoms. Timestep: {1:d}\n".format(Ntr, step)
        g.write( line )
        for i in range( Ntr ):
            line = f.readline()
            step, time, Id, Type, xxx, yyy, zzz, rrr, phi = getTracerData(line)

            #print time
            Type = 1
            if Id in outers:
                Type = 2

            line = "{0:d} {1:g} {2:d} {3:d} {4:g} {5:g} {6:g} {7:g} {8:g} \n".format( step, time, Id, Type, xxx, yyy, zzz, rrr, phi)
            g.write( line )
            print "Writing ", step


def getBinaryPos( time ):

    a     = 1.0
    q     = 1.0
    mu    = q/(1.+q)
    r1   = a*mu
    r2   = a*(1-mu)
    omega = a**( -1.5 )

    phi1_0 = 0.0
    phi2_0 = np.pi
    phi1   = phi1_0 + omega*time
    phi2   = phi2_0 + omega*time

    pl1 = [r1, phi1]
    pl2 = [r2, phi2]

    return pl1, pl2


if __name__ == '__main__':

    stepsize = 25
    cutframe = 800
    infile  = "zbinary2.xyz"
    outfile = "zbin_accret_tot.xyz"
    #find_mini_trs( infile, stepsize, cutframe )
    #tr2mini( infile, outfile )

    reTyping( infile, outfile )
