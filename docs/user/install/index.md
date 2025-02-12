(sec-install)=
# Installation

{numref}`fig:install:choices` provides a guide to select the appropriate method for installing PyLith.
Installation of PyLith on a desktop or laptop machine is, in most cases, very easy.
Binary packages have been created for Linux and macOS platforms.
For Windows users, we recommend installing the Windows Subsystem for Linux and using the Linux binary (see instructions in Section {ref}`sec:install:windows`).
You can also run PyLith inside a Docker container, which provides a virtual Linux environment on any platform that Docker supports, including Linux, macOS, and Windows.
Installation of PyLith on other operating systems - or installation on a cluster - requires building the software from the source code, which can be difficult for inexperienced users.
We have created a small utility called PyLith Installer that makes installing PyLith and all of its dependencies from source much easier.

:::{figure-md} fig:install:choices
<img src="figs/install_choices.*" alt="Guide for selecting the appropriate installation choice based on a hardware and intended use. The installation options are discussed in more detail in the following sections." width = "100%" />

Guide for selecting the appropriate installation choice based on a hardware and intended use.
The installation options are discussed in more detail in the following sections.
:::

Help for installing and using PyLith is available from both a CIG forum and the GitHub issue tracking system <https://github.com/geodynamics/pylith/issues>.
See {ref}`sec-getting-help` for more information.

## Installation of Binary Executable

The binaries are intended for users running on laptops or desktop computers (as opposed to clusters).
The binaries contain the compilers and header files, so users wishing to extend the code can still use the binary and do not need to build PyLith and its dependencies from source.
See {ref}`sec-developer-contributing` for more information on extending PyLith.

Binary executables are available for Linux (glibc 2.12 and later) and macOS (Intel 11.0 and later) from the PyLith web page <https://geodynamics.org/resources/pylith/>.
Users running macOS on MX based computers can use the Intel version.
Users running Windows 10 build 14316 and later can install a Linux bash environment and use the PyLith binary for Linux (see Section {ref}`sec:install:windows` for more information).

:::{tip}
On Linux systems you can check which version of glibc you have by running `ldd-version`

On macOS systems you can check the operating system version by clicking on the Apple icon and *About This Mac*.
:::

### Linux and macOS

1.  Open a terminal window and change to the directory where you want to place the distribution.
    ```{code-block} bash
    $ cd $HOME
    $ mkdir pylith
    $ cd pylith
    ```
2.  Download the Linux or macOS tarball from the PyLith web page
    <https://geodynamics.org/resources/pylith/supportingdocs/>, and save it to
    the desired location, e.g., `$HOME/pylith`.
3.  Unpack the tarball.
    ```{code-block} bash
      # Linux 64-bit
      $ tar -xzf pylith-3.0.1-linux-x86_64.tar.gz

      # macOS
      $ tar -xzf pylith-3.0.1-macOS-11.6.6-x86_64.tar.gz
      ```
4. Set environment variables.
The provided `setup.sh` script only works if you are using bash shell.
If you are using a different shell, you will need to alter how the environment variables are set in `setup.sh`.
```{code-block} bash
$ source setup.sh
Ready to run PyLith.
```

:::{tip}
To bypass macOS quarantine restrictions, simply use command line program `curl` to download the tarball from within a terminal rather than using a web browser.

```
curl -L -O https://github.com/geodynamics/pylith/releases/download/v3.0.1/pylith-3.0.1-macOS-11.6.6-x86_64.tar.gz
```

Alternatively, if you do download the tarball using a web browser, after you unpack the tarball you can remove the macOS quarantine flags using the following commands (requires Administrator access):

```
# Show extended attributes
xattr ./pylith-3.0.1-macOS-11.6.6-x86_64

# Output should be
com.apple.quarantine

# Remove quarantine attributes
sudo xattr -r -d com.apple.quarantine ./pylith-3.0.1-macOS-11.6.6-x86_64
```
:::

:::{warning}
The binary distribution contains PyLith and all of its dependencies.
If you have any of this software already installed on your system, you need to be careful in setting up your environment so that preexisting software does not conflict with the PyLith binary.
By default the `setup.sh` script will prepend to the PATH and PYTHONPATH (for macOS and Linux) and LD_LIBRARY_PATH (for Linux) environment variables.
This will prevent most conflicts.
:::

(sec:install:windows)=
### Windows

PyLith is developed within the Unix/Linux framework, and we do not provide a native PyLith binary distribution for Windows.
The preferred approach to installing PyLith on a computer running Windows is to enable use of a Linux subsystem.
This permits use of the PyLith Linux x86_64 binary within the bash environment.

To enable the Linux subsystem on Windows 10 build 14316 and later (users running an earlier Windows build should use the PyLith Docker container):

1. Go to *Settings* &#8594; *Security* .
2. Under *For developers* select *Developer mode*.
This step should not be required for Windows build 16215 and later.
3. Go to *Control Panel* &#8594; *Programs* &#8594; *Turn Windows Features On or Off*.
4. Enable *Windows Subsystem for Linux* and click *OK* .
5. Restart the computer.
6. Go to *Start* &#8594; *bash*.
You will be prompted to download "Bash on Ubuntu on Windows" from the Windows Store.
Create a user account and password for the bash environment.
7. Install the PyLith Linux x86 binary within the bash environment following the instructions for installing the PyLith binary for Linux.
You will run PyLith within the bash environment just like you would for a Linux operating system.

:::{tip}
To be able to use Gmsh within the Windows Subsystem for Linux, you can install `libglu1-mesa-dev`.
This provides the graphics libraries needed by Gmsh.
:::

### Extending PyLith or Integrating Other Software Into PyLith

:::{note}
New in v3.0.0
:::

We strongly recommend using the [PyLith development environment Docker container](https://pylith-installer.readthedocs.io/en/latest/devenv/index.html) if you want to extend PyLith or integrate PyLith into other software.

## Installation from Source

PyLith depends on a number of other packages (see {numref}`fig:pylith:dependencies`).
This complicates building the software from the source code.
In many cases some of the packages required by PyLith are available as binary packages.
On the one hand, using the binary packages for the dependencies removes the burden of configuring, building, and installing these dependencies, but that can come with its own host of complications if consistent compiler and configuration settings are not used across all of the packages on which PyLith depends.
This is usually not an issue with Linux distributions, such as Fedora, Ubuntu, and Debian that have good quality control; it can be an issue with macOS package managers, such as Fink, MacPorts, and Homebrew, where there is limited enforcement of consistency across packages.
Nevertheless, PyLith can be built on most systems provided the instructions are followed carefully.
PyLith is developed and tested on Linux and macOS.

A small utility, PyLith Installer, removes most of the obstacles in building PyLith and its dependencies from source.
For each package this utility downloads the source code, configures it, builds it, and installs it.
 This insures that the versions of the dependencies are consistent with PyLith and that the proper configure arguments are used.
 The minimum requirements for using the PyLith installer are a C compiler, `tar`, and `wget` or `curl`. Detailed instructions for how to install PyLith using the installer are included in the installer distribution, which is available from the PyLith web page <https://geodynamics.org/resources/pylith/supportingdocs/>.

## Verifying PyLith is Installed Correctly

The easiest way to verify that PyLith has been installed correctly is to run one or more of the examples supplied with the binary and source code.
In the binary distribution, the examples are located in `src/pylith-3.0.1/examples` while in the source distribution, they are located in `pylith-3.0.1/examples`.
{ref}`sec-examples` discusses how to run and visualize the results for the examples.
To run the example discussed in Section {ref}`sec-examples-box-2d`:

```{code-block} bash
$ cd examples/box-2d
$ pylilth step01_axialdisp.cfg
# A bunch of stuff will be written to stdout. The last few lines should be:
 >> .../lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

If you run PyLith in a directory without any input, you will get the error
message:

```{code-block} bash
$ pylith
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.description <- ''
 -- Nonempty string required.
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.arguments <- '[]'
 -- List of command line arguments required.
 >> {default}::
 -- pyre.inventory(error)
 -- metadata.pylith_version <- '[]'
 -- List of PyLith version constraints required.
 >> {default}::
 -- pyre.inventory(error)
 -- meshimporter.meshioascii.filename <- ''
 -- Filename for ASCII input mesh not specified.  To test PyLith, run an example as discussed in the manual.
 >> {default}::
 -- pyre.inventory(error)
 -- timedependent.problem_defaults.name <- ''
 -- Missing required property 'name' in default options for problem.
pylithapp: configuration error(s)
```

This indicates that at a very minimum metadata and the finite-element mesh file must be specified in order to run PyLith.

## Configuration on a Cluster

If you are installing PyLith on a cluster with a batch system, you can configure Pyre such that the `pylith` command automatically submits jobs to the batch queue.
Pyre contains support for the LSF, PBS, SGE, and Globus batch systems.

The command to submit a batch job depends upon the particular batch system used.
Further, the command used in a batch script to launch an MPI program varies from one cluster to the next.
This command can vary between two clusters, even if the clusters use the same batch system!
On some systems, `mpirun` is invoked directly from the batch script.
On others, a special wrapper is used instead.

Properly configured, Pyre can handle job submissions automatically, insulating users from the details of the batch system and the site configuration.

(sec:launchers:schedulers)=
### Launchers and Schedulers

If you have used one of the batch systems, you will know that the batch system requires you to write a script to launch a job.
Fortunately, launching a parallel PyLith job is simplified by Pyre's **launcher** and **scheduler** facilities.
Many properties associated with **launcher** and **scheduler** are pertinent to the cluster you are on, and are best customized in a configuration file.
Your personal PyLith configuration file (`$HOME/.pyre/pylithapp/pylithapp.cfg`) is suitable for this purpose.

Pyre's **scheduler** facility is used to specify the type of batch system you are using (if any):

```{code-block} cfg
[pylithapp]
# The valid values for scheduler are 'lsf", 'pbs', 'globus', and 'none.
scheduler = lsf
# Pyre's launcher facility is used to specify the MPI implementation.
# The valid values for launcher include 'mpich' and 'lam-mpi'.
launcher = mpich
```

You may find the 'dry' option useful while debugging the **launcher** and **scheduler** configuration.
This option causes PyLith to perform a "dry run," dumping the batch script or mpirun command to the console, instead of actually submitting it for execution (the output is only meaningful if you're using a batch system).

```{code-block} bash
# Display the bash script that would be submitted.
$ pylith --scheduler.dry
# Display the mpirun command.
$ pylith --launcher.dry
```

### Running without a Batch System

On a cluster without a batch system, you need to explicitly specify the machines on which the job will run.
Supposing the machines on your cluster are named n001, n002, ..., etc., but you want to run the job on machines n001, n003, n004, and n005 (maybe n002 is down for the moment).
To run an example, create a file named `mymachines.cfg` which specifies the machines to use:

```{code-block} cfg
[pylithapp.launcher]
nodegen = n%03d
nodelist = [1,3-5]
```

The **nodegen** property is a printf-style format string, used in conjunction with **nodelist** to generate the list of machine names.
The `nodelist` property is a comma-separated list of machine names in square brackets.

Now, invoke the following:

```{code-block} bash
$ pylith example.cfg mymachines.cfg
```

This strategy gives you the flexibility to create an assortment of `cfg` files (with one `cfg` file for each machine list) which can be easily paired with different parameter files.

If your machine list does not change often, you may find it more convenient to specify default values for **nodegen** and **nodelist** in `$HOME/.pyre/pylithapp/pylithapp.cfg` (which is read automatically).
Then, you can run any simulation with no additional arguments:

```{code-block} console
$ pylilth example.cfg
```

:::{warning}
This assumes your machine list has enough nodes for the simulation in question.
:::

You will notice that a machine file `mpirun.nodes` is generated.
It will contain a list of the nodes where PyLith has run.

### Using a Batch System

Many clusters use some implementation of a PBS (e.g., TORQUE/Maui) or LSF batch system.
The examples below illustrate use of some of the more important settings.
You may need to make use of more options or adjust these to submit jobs on various cluster.
These settings are usually placed in `$HOME/.pyre/pylithapp/pylithapp.cfg` or in a system-wide configuration file.
They can be overridden on the command line, where one typically specifies the number of compute nodes and number of processes per compute node, the job name, and the allotted time for the job:

```{code-block} bash
$ pylith example1.cfg \
--job.queue=debug \
--job.name=example1 \
--job.stdout=example1.log \
--job.stderr=example1.err \
--job.walltime=5\*minute \
--nodes=4
```

:::{important}
The value for nodes is equal to the number of compute nodes times the number of processes (usually the number of cores) requested per compute node.
Specifying the number of processes per compute node depends on the batch system.
For more information on configuring Pyre for your batch system, see CIG's Pythia page <https://geodynamics.org/resources/pythiaXX/>.
:::

#### LSF Batch System

```{code-block} cfg
[pylithapp]
scheduler = lsf
# the type of batch system

[pylithapp.lsf]
bsub-options = [-a mpich_gm]
# special options for 'bsub'

[pylithapp.launcher]
command = mpirun.lsf
'mpirun' command to use on our cluster

[pylithapp.job]
queue = normal
# default queue for jobs
```

#### PBS Batch System

```{code-block} cfg
[pylithapp]
scheduler = pbs
# the type of batch system

[pylithapp.pbs]
shell = /bin/bash
# submit the job using a bash shell script

# Export all environment variables to the batch job
# Send email to johndoe@mydomain.org when the job begins, ends, or aborts
qsub-options = -V -m bea -M johndoe@mydomain.org

[pylithapp.launcher]
command = mpirun -np ${nodes} -machinefile ${PBS_NODEFILE}
```

For most PBS batch systems you can specify N processes per compute node via the command line argument `--scheduler.ppn=N`.

% End of file