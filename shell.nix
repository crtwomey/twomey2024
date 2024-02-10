# nix development environment with project dependencies

{ pkgs ? import <nixpkgs> {} }:

with pkgs;

let
  performance = rPackages.buildRPackage {
    name = "performance";
    # need latest version to fix problem with lmer models.
    src  = fetchFromGitHub {
      owner  = "easystats";
      repo   = "performance";
      rev    = "90956a9dfacd0d19a12d4a618bdbbc1a7ac665b0";
      sha256 = "0lw1q0mw3v5phfascl5lwzd8myvacjdzp9gi1pfm2w0q1xph42zw";
    };
    buildInputs = [
      R
      rPackages.bayestestR
      rPackages.insight
    ];
  };
  emdist = rPackages.buildRPackage {
    name = "emdist";
    # need version 3-2 for setting the maximum number of iterations, but this
    # version is not yet available on CRAN. Pull from github instead.
    src  = fetchFromGitHub {
      owner  = "s-u";
      repo   = "emdist";
      rev    = "cc68101a6606aa7f17ce43ad0cea13071c7ba7f0";
      sha256 = "0aiqr7pn95xwk7bg4nxbkd22jx2m8pin5amnkg7319f2lv18pgyk";
    };
  };
  ggsankey = rPackages.buildRPackage {
    name = "ggsankey";
    # need latest version to fix problem with lmer models.
    src  = fetchFromGitHub {
      owner  = "davidsjoberg";
      repo   = "ggsankey";
      rev    = "be08dd0f86eaee9f9ff9e7ff95d47930660a3c36";
      sha256 = "0acpmydqqc91pq5p9wpkpmgqp3nhiljabd7d3i00kwhjxgm2bvba";
    };
    buildInputs = with rPackages; [
      R ggplot2 dplyr stringr purrr tidyr magrittr
    ];
  };
in mkShell {
  packages = [
    libintl

    # Markdown
    pandoc # generate pdf from markdown

    # R environment
    R
    rPackages.targets        # reproducible builds
    rPackages.tarchetypes    # knitr target
    rPackages.visNetwork     # visualize targets
    rPackages.pbapply        # progress bar
    rPackages.parallelly     # parallel execution

    rPackages.tidyverse      # tidy data wrangling
    rPackages.paletteer      # color palettes
    rPackages.viridis        # accessible color gradients
    rPackages.knitr          # RMarkdown
    rPackages.tikzDevice     # use tikz from R

    rPackages.colorscience   # color space conversions
    rPackages.lpSolve        # linear assignment solveq
    rPackages.igraph         # graph modularity 
    rPackages.deldir         # voronoi tesselation
    rPackages.geometry       # n-dimensional convex hull
    emdist                   # earth mover's distance
    rPackages.gtools         # dirichlet distribution
    rPackages.ks             # kernel density estimation

    rPackages.lme4           # mixed effects models
    rPackages.DHARMa         # model diagnostics
    performance              # model performance
    rPackages.see            # performance dependency
    rPackages.ggrepel        # performance dependency
    rPackages.qqplotr        # performance dependency
    rPackages.ggeffects      # partial residuals

    rPackages.sjPlot         # plot model
    rPackages.cowplot        # plot composition
    rPackages.patchwork      # plot composition
    rPackages.ape            # dendrogram plots
    rPackages.plot3D         # 3D plotting
    rPackages.quantreg       # geom_boxplot with weights
    rPackages.ggtext         # images as labels
    rPackages.figpatch       # images in patchwork
    ggsankey                 # sankey diagram
    rPackages.pdftools       # add ext pdf figures
  ];
}
