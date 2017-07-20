import sys
import discopy as dp

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("usage: $ python remapCheckpoint.py t <checkpoint.h5 ...>")
        print("Sets time in checkpoints to t.")
        sys.exit()

    t = float(sys.argv[1])

    for f in sys.argv[2:]:
        print("Processing {0:s}".format(f))
        dp.setTime(f, t)
