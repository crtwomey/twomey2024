[![DOI](https://zenodo.org/badge/755312042.svg)](https://zenodo.org/doi/10.5281/zenodo.10642708)

# Code and Data Supplementary Information for Twomey et al. (2024)

This repository contains supplementary code and data for Twomey et al. (2024).

Project dependencies are listed in the `shell.nix` file. An environment
containing these dependencies can be constructed by using the `nix` package
manager and running the `nix-shell` command from the top-level directory.

The `_targets.R` file specifies a
[targets](https://github.com/ropensci/targets) pipeline. Launching an `R` REPL
and running the following from the top-level directory will reproduce the
results for Twomey et al. 2024:

```R
library(targets)
tar_make()
```

The data in this repository are included for reproducibility. They are derived
from the [WCS Data Archives](https://www1.icsi.berkeley.edu/wcs/data.html) and
[Twomey et al. 2021](https://doi.org/10.1073/pnas.2109237118). If you use this
data in another work, please ensure the original source is cited.
