from distutils.core import setup, Extension
import numpy.distutils.misc_util


calc_dir = "../../Calc"
inc = numpy.distutils.misc_util.get_numpy_include_dirs()
inc.append(calc_dir)

setup(
    ext_modules=[Extension("_calc", ["_calc.c",calc_dir+"/bondi.c"])],
    include_dirs=inc
)
