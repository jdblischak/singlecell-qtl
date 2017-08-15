# singlecell-qtl

A [workflowr][] project.

## Setup

To ensure all contributors are using the same computational environment, we use 
[conda][] to manage software dependencies (made possible by the [bioconda][] and
[conda-forge][] projects). Please complete the following steps to replicate the
computing environment. Note that this is only guaranteed to work on a Linux-64
based architecture, but in theory should be able to work on macOS as well.

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

1. Download and install the R package [workflowr][] from source (unfortunately
you can't use `devtools::install_github()` because of incompatibilities
introduced in the conda environment)
    ```
    wget -O /tmp/workflowr-v0.7.0.tar.gz https://github.com/jdblischak/workflowr/archive/v0.7.0.tar.gz
    R CMD INSTALL /tmp/workflowr-v0.7.0.tar.gz
    ```

If there are updates to `environment.yaml`, you can update the "scqtl"
environment by running `conda env udpate --file environment.yaml`.

[bioconda]: https://bioconda.github.io
[conda]: https://conda.io/docs/
[conda-forge]: https://conda-forge.org/
[workflowr]: https://github.com/jdblischak/workflowr
