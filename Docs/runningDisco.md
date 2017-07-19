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

### Running Disco ###

Disco reads in runtime parameters from a file called `in.par` which must exist in the working directory where you run the executable.  A generic template `in.par` exists in the source directory as `in.par.template`. Copy this into `in.par` and change as you wish.  

#### Using Pre-built Templates ####

We have provided several template runs (consisting of a `Makefile_opt.in` and an `in.par`) in the `Templates/` directory.  Copy them yourself, or run `make <template name>` to quickly clean the directory, copy the template files, and compile Disco.  For example to run the vortex test simply type:

```bash
$ make vortex
$ ./disco
```

Most Disco features are utilized in at least one of the templates. These provide a good starting point to customize your own runs of Disco.
