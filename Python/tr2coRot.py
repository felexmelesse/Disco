import numpy as np

def tr2coRot( inFile, outFile ):

    f = open( inFile , "r" )
    g = open( outFile, "w" )

    a  = 1.0
    om = pow( a, -1.5 )
    step = 0

    line = "start"
    while line != " ":
        line = f.readline()
        Ntr  = int( line )
        skip = f.readline()
        line = "{0:d} \nAtoms. Timestep: {1:d}\n".format(Ntr, step)
        g.write( line )
        for i in range( Ntr ):
            line = f.readline()
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

            print phi, om, time
            phi = phi - om*time
            xxx = rrr*np.cos(phi)
            yyy = rrr*np.sin(phi)

            line = "{0:g} {1:d} {2:d} {3:d} {4:g} {5:g} {6:g} {7:g} {8:g} \n".format(time, step, Id, Type, xxx, yyy, zzz, rrr, phi)
            g.write( line )
            print step


if __name__ == '__main__':

    infile  = "zbinary.xyz"
    outfile = "zbin_corot.xyz"
    tr2coRot( infile, outfile )
