# singlecell-qtl

A [workflowr][] project.

## Setup

To ensure all contributors are using the same computational environment, we use
[conda][] to manage software dependencies (made possible by the [bioconda][] and
[conda-forge][] projects). Please complete the following steps to replicate the
computing environment. Note that this is only guaranteed to work on a Linux-64
based architecture, but in theory should be able to work on macOS as well. All
commands shown below are intended to be run in a Bash shell from the root of the
project directory.

1. Install Git and register for an account on GitHub

1. Download and install Miniconda ([instructions](https://conda.io/miniconda.html))

1. Clone this repository (or your personal fork) using `git clone`

1. Create the conda environment "scqtl" using `environment.yaml`
    ```
    conda env create --file environment.yaml
    ```

1. To use the conda enviroment, you must first activate it by running `source
activate scqtl`. This will override your default settings for R, Python, and
various other software packages. When you are done working on this project, you
can either logout of the current session or deactivate the environment by
running `source deactivate`.

1. Initialize git-lfs and download latest version of large data files
    ```
    git lfs install
    git lfs pull
    ```

If there are updates to `environment.yaml`, you can update the "scqtl"
environment by running `conda env udpate --file environment.yaml`.

**Warning:** If you are using RStudio, you need to ensure that it recognizes
  your conda environment. If you launch RStudio by clicking on an icon, it
  doesn't use the current environment you have configured in your shell. On a
  Linux-based system, the solution is to launch RStudio directly from the shell
  with `rstudio`. On macOS, running `open -a rstudio` from the shell causes
  RStudio to recognize most of the environemnt variables, but strangely it does
  not set the correct library path to the conda R packages. Suggestions for how
  to fix this are welcome.

[bioconda]: https://bioconda.github.io
[conda]: https://conda.io/docs/
[conda-forge]: https://conda-forge.org/
[workflowr]: https://github.com/jdblischak/workflowr
