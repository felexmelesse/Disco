# Running Disco #

### First Time Setup ###

Disco relies on the HDF5 (for efficient i/o) and MPI (for parallism) libraries. The location of these libraries, as well as other system-specifc options, are specified by `Makefile_dir.in`. We provide a template in `Makefile_dir.in.template`. First copy the template into a new file:

```bash
$ cp Makefile_dir.in.template Makefile_dir.in
```

Then modify the entries as necessary. You should not need to change this file much (if at all) once you have Disco running.

To compile a Disco binary you also need a `Makefile_opt.in` which sets run-specific options like the initial condition, boundary conditions, and system of equations to solve.  A template exists in `Makefile_opt.in.template`, copy it into a new file and modify as you need:

```bash
$ cp Makefile_opt.in.template Makefile_opt.in
```

With these you should be able to compile Disco by just running `make`:

```bash
$ make
```

Voila! You have a Disco binary to run.

#### Setting up *discopy* ####

Disco includes a Python package called `discopy` to interface with the data files, make geometry-aware plots, and perform simple analysis. The package is compatible with both Python2 and Python3.

To install `discopy` navigate to the `Python` directory and run the `setup.py` file.  This will install `discopy` into your system `site-packages` and allow you to import `discopy` from anywhere on your system.  

```bash
$ cd Python/
$ python setup.py
```

If you plan on modifying `discopy` consider adding the `develop` argument to your `setup.py` call, this only installs symlinks into the system `site-packages` directory, allowing your changes to instantly be reflected in the module without having to re-run `setup.py`.

```bash
$ python setup.py develop
```

### Running Disco ###

Disco reads in runtime parameters from a file called `in.par` which must exist in the working directory where you run the executable.  A generic template `in.par` exists in the source directory as `in.par.template`. Copy this into `in.par` and change as you wish.  

#### Using Pre-built Templates ####

We have provided several template runs (consisting of a `Makefile_opt.in` and an `in.par`) in the `Templates/` directory.  Copy them yourself, or run `make <template name>` to quickly clean the directory, copy the template files, and compile Disco.  For example to run the vortex test simply type:

```bash
$ make vortex
$ ./disco
```

Most Disco features are utilized in at least one of the templates. These provide a good starting point to customize your own runs of Disco.
