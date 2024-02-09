#
# _targets.R
#
# Pipeline for reproducing Twomey et al. 2024.
# All code is open source under the GPL v3 (see LICENSE)
#

library(targets)
library(tarchetypes)

# enable progress bars for long computations
library(pbapply)
pboptions(type = "txt", char = "-")

# whether or not to generate the manuscript
skip_manuscript_targets <- FALSE

# use tikz for plotting
tikz_enabled <- TRUE

# parameters and utilities
source("src/common/definitions.R")
source("src/common/plotting.R")
source("src/common/parallel.R")

# load WCS data
source("src/wcs_data.R")
source("src/wcs_plot.R")
source("src/wcs_distance.R")

# word clusters
source("src/word_cluster_generation.R")
source("src/word_cluster_assignments.R")
source("src/word_cluster_plot.R")

# rate-distortion
source("src/rd.R")
source("src/rd_example.R")

# n-1+1 vocabularies
source("src/vocabulary_generation.R")
source("src/vocabulary_analysis.R")
source("src/vocabulary_summary_plot.R")
source("src/vocabulary_examples_plot.R")
source("src/vocabulary_likelihoods_check.R")

# n to n+1 transitions
source("src/word_transitions.R")
source("src/word_transitions_plot.R")

# paths
data_dir        <- "data"        # data inputs
fig_dir         <- "fig"         # figure outputs
manuscript_dir  <- "manuscript"  # manuscript text
publication_dir <- "publication" # publication outputs

# global packages
tar_option_set(packages = c(
  "tidyverse",
  "colorscience",
  "parallelly",
  "listenv",
  "emdist",
  "lpSolve",
  "igraph",
  "patchwork",
  "paletteer",
  "viridis",
  "ggsankey",
  "ggtext",
  "figpatch"
))

#------------------------------------------------------------------------------
# data targets
#
# defines a list of data aggregation targets for the World Color Survey. The
# targets load and aggregates data from multiple files, including color chip
# references, vocabulary data, language identifiers, and inferences. The final
# targets are the aggregated data set and color terms built from every language
# in the survey.
#
data_aggregation_targets <- list(
  # WCS color chip reference
  tar_target(wcs_chips_file,
    file.path(data_dir, "wcs-chips.rds"),
    format = "file"
  ),

  # pre-processed WCS vocabulary data
  tar_target(wcs_vocabularies_file,
    file.path(data_dir, "wcs-vocabularies.rds"),
    format = "file"
  ),

  # (corrected) WCS language identifiers
  tar_target(wcs_info_file,
    file.path(data_dir, "WCS_SIL_codes.txt"),
    format = "file"
  ),

  # inferences for WCS data from Twomey et al. 2021
  tar_target(wcs_infs_file,
    file.path(data_dir, "wcs-inf.rds"),
    format = "file"
  ),

  # load and aggregate data
  tar_target(wcs, load_wcs_data(
    wcs_chips_file,
    wcs_vocabularies_file,
    wcs_info_file,
    wcs_infs_file
  )),

  # load color terms from every WCS language
  tar_target(wcs_terms, build_wcs_terms(wcs))
)

#------------------------------------------------------------------------------
# word cluster targets
#
# defines a list of targets for building word clusters for the World Color
# Survey. It includes targets for computing distances and adjacency matrices
# between terms, generating clusters, mapping terms to words, computing
# vocabulary size frequencies, and generating a function for computing word
# assignments.
#
term_clustering_targets <- list(
  # generate function for computing earth mover's distance for distributions
  # over the WCS stimulus set.
  tar_target(wcs_distance, build_wcs_distance_function(wcs)),

  # compute pairwise distances between all wcs terms
  tar_target(term_distances, compute_term_distances(
    wcs_terms,
    wcs_distance
  )),

  # compute pairwise adjacency matrix between all wcs terms
  tar_target(term_adjacency, compute_term_adjacency(
    wcs_terms,
    term_distances
  )),

  # compute clusters based on pairwise adjacencies
  tar_target(word_clusters, compute_word_clusters(
    wcs,
    wcs_terms,
    term_adjacency
  )),

  # words for each language (mapping of terms to words and removal of synonyms)
  tar_target(wcs_words, compute_wcs_words(
    wcs,
    wcs_terms,
    word_clusters
  )),

  # frequency of each vocabulary size across the WCS languages
  tar_target(wcs_nword_frequencies, compute_nword_frequencies(
    wcs_words
  )),

  # generate function for computing word assignments
  tar_target(assign_word, build_word_assignment_function(
    wcs_distance,
    word_clusters
  ))
)

#------------------------------------------------------------------------------
# targets for generating n+1, n-1, and n-1+1 results
#
# defines targets for generating different types of vocabularies and
# calculating various statistics related to these vocabularies. These targets
# include generating rate-distortion optimal vocabularies, generating n+1 and
# n-1 word vocabularies, generating de novo vocabularies, computing pairwise
# distances between vocabularies, approximating the likelihood of an n-1
# vocabulary being the ancestor of a vocabulary, and computing the absolute and
# effective number of modes and RD efficiencies for different types of
# vocabularies.
#
word_evolution_targets <- list(
  # rate-distortion optimal vocabularies
  tar_target(rd_optimal_vocabularies, rd_optimize_vocabularies(
    wcs, wcs_words
  )),

  # for each language, generate a set of rate-distortion optimal n+1 word
  # vocabularies.
  tar_target(nplus1_vocabularies, generate_nplus1_rd_vocabularies(
    wcs, rd_optimal_vocabularies
  )),

  # for each language, generate a set of rate-distortion optimal n-1 word
  # vocabularies corresponding to the original n word vocabulary with one
  # word deleted.
  tar_target(nminus1_vocabularies, generate_nminus1_vocabularies(
    rd_optimal_vocabularies
  )),

  # for each potential n-1 word vocabulary for each language, generate a set
  # of potential n word vocabularies (n-1+1 vocabularies).
  tar_target(nminus1plus1_vocabularies, generate_nplus1_vocabularies(
    wcs, nminus1_vocabularies
  )),

  # for each language, generate a set of de novo vocabularies of the same
  # number of words and using the same distribution of communicative needs,
  # p(x), and beta.
  tar_target(denovo_vocabularies, generate_denovo_vocabularies(
    wcs, wcs_words
  )),

  # for each language, for each de novo mode for that language, generate a set
  # of potential n+1 word languages.
  tar_target(denovo_nplus1_vocabularies, generate_nplus1_vocabularies(
    wcs, denovo_vocabularies
  )),

  # quantification of variability in probabilities of successor vocabularies
  tar_target(nplus1_variability, compute_nplus1_variability(
    nplus1_vocabularies,
    denovo_nplus1_vocabularies
  )),

  # generate de novo n-1 vocabularies with the same communicative needs and
  # beta as the original size n vocabulary, but with a size of n-1 and no
  # initial positioning based on the original size n vocabulary.
  tar_target(denovo_nminus1_vocabularies, generate_denovo_vocabularies(
    wcs, wcs_words, -1
  )),

  # n+1 de novo vocabularies (note: different from denovo_nplus1_vocabularies.
  # Above we make a denovo n vocabulary then add a term. Here we make a denovo
  # n+1 vocabulary directly.)
  tar_target(nplus1_denovo_vocabularies, generate_denovo_vocabularies(
    wcs, wcs_words, 1
  )),

  # compute distances between all pairs of n+1 vocabularies
  # (historical and de novo).
  tar_target(nplus1_pairwise_distances, compute_nplus1_pairwise_distances(
    nplus1_vocabularies,
    nplus1_denovo_vocabularies,
    wcs_distance
  )),

  # compute pairwise distances between n+1 vocabularies based on the minimum
  # Wasserstein distance assignment between them.
  tar_target(nplus1_min_distances, compute_nplus1_min_distances(
    nplus1_pairwise_distances
  )),

  # generate n-1+1 de novo vocabularies (de novo n-1 vocabularies with 1 word
  # subsequently added).
  tar_target(denovo_nminus1plus1_vocabularies, generate_nplus1_vocabularies(
    wcs, denovo_nminus1_vocabularies
  )),

  # n-1+1 distances to original vocabulary
  tar_target(denovo_nminus1plus1_distances,
    compute_denovo_nminus1plus1_distances(
      wcs_words,
      denovo_nminus1plus1_vocabularies,
      wcs_distance
    )
  ),

  # n-1+1 pairwise distances to other de novo modes
  tar_target(denovo_nminus1plus1_pairwise_distances,
    compute_denovo_nminus1plus1_pairwise_distances(
      denovo_nminus1plus1_vocabularies,
      wcs_distance
    )
  ),

  # table of distances and frequencies for de novo n-1+1 vocabularies
  tar_target(denovo_nminus1plus1_fd_table,
    build_frequency_and_distance_table(
      denovo_nminus1plus1_distances,
      denovo_nminus1plus1_vocabularies
    )
  ),

  # determine reasonable bounds for the prior over the scale parameter
  tar_target(lb_distance_scales,
    compute_scale_parameter_prior_bounds(denovo_nminus1plus1_fd_table)
  ),

  # approximate the likelihood of an n-1 vocabulary being the ancestor of a
  # vocabulary based on the distance of +1 versions of the n-1 vocabulary to
  # the focal n word vocabulary.
  tar_target(denovo_nminus1plus1_likelihoods,
    compute_precursor_likelihoods(
      denovo_nminus1plus1_fd_table,
      lb_distance_scales
    )
  ),

  # compute the absolute and effective number of modes for both historical and
  # de novo n+1 vocabularies.
  tar_target(modes_summary, summarize_vocabulary_modes(
    nplus1_vocabularies,
    nplus1_denovo_vocabularies
  )),

  # compute RD efficiencies for n+1 historical and de novo vocabulary modes
  tar_target(modes_efficiency, compute_vocabulary_modes_efficiency(
    nplus1_vocabularies,
    nplus1_denovo_vocabularies
  ))
)

#------------------------------------------------------------------------------
# targets for summarizing vocabulary modes
#
# builds tables of mode numbers and frequencies for historical and de novo n
# and n+1 vocabularies. Then computes color word assignments for every term in
# every mode. E.g.
#
# > head(word_modes_table)
# A tibble: 6 × 13
#   language type       parent_mode_id mode_id mode_nterms mode_freq mode
#      <int> <chr>               <dbl>   <dbl>       <int>     <dbl>
# 1        1 historica…             NA       1           6         1
# 2        1 historica…             NA       1           6         1
# 3        1 historica…             NA       1           6         1
# 4        1 historica…             NA       1           6         1
# 5        1 historica…             NA       1           6         1
# 6        1 historica…             NA       1           6         1
# ... with 4 more variables: px_xhat <list>, avg_color <chr>, word <int>
#
# if the vocabulary mode was generated from another vocabulary mode,
# parent_mode_id records the id of the parent (e.g. a n+1 word vocabulary
# generated from an n word vocabulary). The word column contains the assignment
# of the term to color word category.
#
modes_table_targets <- list(
  # helpder function for building term mode tables
  tar_target(table_builder, term_modes_table_builder(wcs)),

  # historical n table
  tar_target(historical_n_term_modes_table,
    build_historical_n_term_modes_table(
      table_builder, rd_optimal_vocabularies
    )
  ),

  # historical n+1 table
  tar_target(historical_nplus1_term_modes_table,
    build_historical_nplus1_term_modes_table(
      table_builder, nplus1_vocabularies
    )
  ),

  # de novo n table
  tar_target(denovo_n_term_modes_table,
    build_denovo_n_term_modes_table(
      table_builder, denovo_vocabularies
    )
  ),

  # de novo n+1 table
  tar_target(denovo_nplus1_term_modes_table,
    build_denovo_nplus1_term_modes_table(
      table_builder, denovo_nplus1_vocabularies
    )
  ),

  # vocabulary modes table summarizing number of modes and their frequency for
  # the four (historical n, n+1, and de novo n, n+1) tables constructed above.
  tar_target(term_modes_table, build_term_modes_table(
    historical_n_term_modes_table,
    historical_nplus1_term_modes_table,
    denovo_n_term_modes_table,
    denovo_nplus1_term_modes_table
  )),

  # compute assignments of terms to words
  tar_target(word_modes_table, compute_modes_table_word_assignments(
    assign_word,
    term_modes_table
  ))
)

#------------------------------------------------------------------------------
# targets for n to n+1 transition results
#
# defines targets for summarizing n to n+1 word transition results, including
# building a table of connections between parent and child terms, reformatting
# the n->n+1 table as a table of transitions, reformatting the transition table
# as a set of transition matrices, and computing conditional transition
# probabilities.
#
word_transition_targets <- list(
  # builds a table of connections between parent and child terms
  tar_target(n_to_nplus1_table, build_n_to_nplus1_table(
    word_modes_table
  )),

  # reformats the n -> n+1 table as a table of transitions
  tar_target(n_to_nplus1_transitions, build_n_to_nplus1_transitions(
    n_to_nplus1_table
  )),

  # reformats the transition table as a set of transition matricies
  tar_target(joint_transition_matricies, build_joint_transition_matricies(
    n_to_nplus1_transitions
  )),

  # computes conditional transition probabilities
  tar_target(conditional_transition_matricies,
    build_conditional_transition_matricies(
      n_to_nplus1_transitions
    )
  )
)

#------------------------------------------------------------------------------
# fig 1 - WCS color category dictionary
#
fig_1_targets <- list(
  # reference plot of the WCS color stimuli
  tar_target(wcs_color_stimuli_plot, plot_wcs_color_stimuli(
    wcs, point_size = 3, axis_text_size = 6
  )),

  # show the mappings for each word (cluster of terms)
  tar_target(word_plots, plot_words(
    wcs,
    word_clusters
  )),

  # histogram of vocabulary size for WCS languages
  tar_target(distribution_of_nwords_figure, plot_distribution_of_nwords(
    wcs_nword_frequencies
  )),

  # combine panels for full figure
  tar_target(word_cluster_figure, plot_word_cluster_figure(
    wcs_color_stimuli_plot,
    distribution_of_nwords_figure,
    word_plots
  )),

  tar_target(word_cluster_fig, fig(
      word_cluster_figure, "word_clusters", 8, 7
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig 2 - evidence for semantic shifts based on WCS
#
fig_2_targets <- list(
  # show conditional transition panels for transitions between the first four
  # stages (3->4, 4->5, 5->6, and 6->7)
  tar_target(conditional_transitions_figure, plot_conditional_transitions(
    word_clusters,
    conditional_transition_matricies
  )),

  # sankey diagram of transitions between vocabulary size stages
  tar_target(sankey_transitions_figure, plot_sankey_transitions(
    word_clusters,
    n_to_nplus1_transitions
  )),

  # show transition matricies and sankey diagram side-by-side
  tar_target(transitions_figure, plot_transitions_figure(
    conditional_transitions_figure,
    sankey_transitions_figure
  )),

  tar_target(transitions_fig, fig(
      transitions_figure, "transitions", 10, 9
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig 3 - example vocabulary successor modes
#
fig_3_targets <- list(
  # show an example comparison between a WCS language's distribution over
  # likely n+1 vocabularies and a related de novo n+1 distribution
  tar_target(historical_denovo_comparison_figure,
    plot_historical_denovo_comparison(
      wcs,
      rd_optimal_vocabularies,
      denovo_vocabularies,
      nplus1_vocabularies,
      denovo_nplus1_vocabularies,
      word_clusters,
      n_to_nplus1_transitions
    )
  ),

  tar_target(comparison_fig, fig(
      historical_denovo_comparison_figure, "comparison", 9, 4.75
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig 4 - succesor modes as a function of vocabulary size
#
fig_4_targets <- list(
  # number of n+1 vocabulary modes as a function of vocabulary size
  tar_target(nmodes_summary_figure, plot_nmodes_summary(
    modes_summary
  )),

  # effective number of n+1 vocabulary modes as a function of vocabulary size
  tar_target(eff_nmodes_summary_figure, plot_eff_nmodes_summary(
    modes_summary
  )),

  # rate-distortion tradeoff across all n+1 word vocabulary modes
  tar_target(modes_efficiency_figure, plot_modes_efficiency(
    modes_efficiency
  )),

  tar_target(postprocessed_nplus1_pairwise_distances,
    compute_postprocessed_nplus1_pairwise_distances(
      wcs_words,
      nplus1_min_distances,
      nplus1_pairwise_distances
    )
  ),

  tar_target(pairwise_distance_percent_difference,
    compute_pairwise_distance_percent_difference(
      postprocessed_nplus1_pairwise_distances
    )
  ),

  # pairwise minimum distance between n+1 word modes
  tar_target(nplus1_min_distances_figure, plot_nplus1_min_distances(
    postprocessed_nplus1_pairwise_distances
  )),

  # combined plot of vocabulary mode summaries
  tar_target(modes_summary_figure, plot_modes_summary(
    nmodes_summary_figure,
    eff_nmodes_summary_figure,
    modes_efficiency_figure,
    nplus1_min_distances_figure
  )),

  tar_target(modes_summary_fig, fig(
      modes_summary_figure, "modes_summary", 10, 10
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig 5 - ancestral language reconstructions
#
fig_5_targets <- list(
  # show examples of inferred ancestral vocabularies
  tar_target(nminus1_vocabularies_figure, plot_nminus1_vocabularies(
    wcs, wcs_words,
    rd_optimal_vocabularies,
    denovo_nminus1_vocabularies,
    denovo_nminus1plus1_likelihoods
  )),

  tar_target(nminus1_vocabularies_fig, fig(
      nminus1_vocabularies_figure, "nminus1_vocabularies", 12.5, 4
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig A1 - RD example
#
fig_a1_targets <- list(
  # generate RDFC results for a simple example with 9 points in a grid
  tar_target(rd_example, gen_example_rd_tradeoff()),

  # plot the example rate-distortiont tradeoff curves
  tar_target(rd_example_figure, plot_rd_example(rd_example)),

  tar_target(rd_example_fig, fig(
      rd_example_figure, "rd_example", 5, 4, use_tikz = tikz_enabled
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# fig A2 - sensitivity check
#
fig_a2_targets <- list(
  # sensitivity of likelihoods to choice of scale parameter
  tar_target(likelihood_sensitivity, compute_likelihood_sensitivity(
    denovo_nminus1plus1_fd_table,
    denovo_nminus1plus1_likelihoods
  )),

  # example of varying scale parameter
  tar_target(sensitivity_example, compute_sensitivity_example(
    denovo_nminus1plus1_fd_table
  )),

  # combined sensitivity plots
  tar_target(sensitivity_figure, plot_likelihood_sensitivities(
    likelihood_sensitivity,
    sensitivity_example,
    lb_distance_scales
  )),

  tar_target(sensitivity_fig, fig(
      sensitivity_figure, "sensitivity", 10, 3, use_tikz = FALSE
    ), format = "file"
  )
)

#------------------------------------------------------------------------------
# manuscript
#
manuscript_targets <- list(
  # package together results for use in the manuscript
  tar_target(results, {
      tar_load(word_clusters)
      tar_load(wcs_nword_frequencies)
      tar_load(nplus1_variability)
      tar_load(denovo_nminus1plus1_likelihoods)
      tar_load(pairwise_distance_percent_difference)
      example_likelihoods <- lapply(fig5_language_examples, function(lnum) {
        get_language_likelihoods(lnum, denovo_nminus1plus1_likelihoods)
      })
      output_file <- file.path(manuscript_dir, "results.rds")
      saveRDS(list(
          nword_clusters                = length(word_clusters$grp),
          wcs_nword_frequencies         = wcs_nword_frequencies,
          initial_frequency_of_new_term = initial_frequency_of_new_term,
          term_adjacency_distance_scale = term_adjacency_distance_scale,
          num_gen_denovo_modes          = num_gen_denovo_modes,
          same_mode_threshold           = same_mode_threshold,
          nplus1_variability            = nplus1_variability,
          likelihood_change_threshold   = likelihood_change_threshold,
          example_likelihoods           = example_likelihoods,
          pairwise_distances            = pairwise_distance_percent_difference
        ), file = output_file
      )
      output_file
    },
    format = "file"
  ),

  # generate tex file
  tar_knit(manuscript, {
      tar_cancel(skip_manuscript_targets)
      file.path(manuscript_dir, "main.Rtex")
    },
    output = file.path(publication_dir, "main.tex")
  ),

  # copy definitions tex file to publication directory
  tar_force(tex_definitions, {
      tar_cancel(skip_manuscript_targets)
      definitions_file <- file.path(manuscript_dir, "definitions.tex")
      output_file      <- file.path(publication_dir, "definitions.tex")
      system(paste("cp", definitions_file, output_file))
      output_file
    },
    format = "file",
    force  = TRUE
  ),

  # copy bib tex file to publication directory
  tar_force(tex_bib, {
      tar_cancel(skip_manuscript_targets)
      bib_file    <- file.path(manuscript_dir, "refs.bib")
      output_file <- file.path(publication_dir, "refs.bib")
      system(paste("cp", bib_file, output_file))
      output_file
    },
    format = "file",
    force  = TRUE
  ),

  # copy bib sty file to publication directory
  tar_force(tex_bst, {
      tar_cancel(skip_manuscript_targets)
      bst_file    <- file.path(manuscript_dir, "refs.bst")
      output_file <- file.path(publication_dir, "refs.bst")
      system(paste("cp", bst_file, output_file))
      output_file
    },
    format = "file",
    force  = TRUE
  ),

  # generate pdf from tex file
  tar_force(publication, {
      tar_cancel(skip_manuscript_targets)
      tar_load(tex_definitions)
      tar_load(tex_bib)
      tar_load(tex_bst)
      tar_load(results)
      tar_load(word_cluster_fig)
      tar_load(transitions_fig)
      tar_load(comparison_fig)
      tar_load(modes_summary_fig)
      tar_load(nminus1_vocabularies_fig)
      tar_load(rd_example_fig)
      tar_load(sensitivity_fig)
      tar_load(manuscript)
      system(
        paste0("latexmk -bibtex -pdf",
          " -output-directory=\"", publication_dir, "\"",
          " -pdflatex=\"xelatex %O",
            " -file-line-error",
            " -interaction=nonstopmode",
            " -synctex=1 %S",
          "\" ",
          manuscript[1]
        )
      )
      paste0(tools::file_path_sans_ext(manuscript[1]), ".pdf")
    },
    format = "file",
    force  = TRUE
  )
)

#------------------------------------------------------------------------------
# all targets
#
c(
  data_aggregation_targets,
  term_clustering_targets,
  word_evolution_targets,
  modes_table_targets,
  word_transition_targets,
  fig_1_targets,
  fig_2_targets,
  fig_3_targets,
  fig_4_targets,
  fig_5_targets,
  fig_a1_targets,
  fig_a2_targets,
  manuscript_targets
)
