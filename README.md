# 2019_dsb_NHEJ

Brief description of files included in this project:

environment_functions.R
This contains all the functions, libraries and imports external data sets used for analyzing and plotting data. Run this code before attempting to run any of the other plotting code.

MNase_coverage_plot.R
This code recreates the MNase coverage plots for figure 1 in the paper.

idealized_nuc_profiling.R
This code is used to profile what an idealized nucleosome looks like in a given MNase-seq experiment. It will generate parameters for an idealized kernel that can be used for the subsequent cross correlation analysis and dot plots.

dot_plots_2d_cross_correlation_pictographs.R
This code will plot both the dot plots (along with 2d cross correlation trace) and pictographs summarize the chromatin dot plot. Be sure to adjust the start and end coordinates if you want to look at another locus and to point this code to use sacCer3 files instead of chr2 files if you want to look at other chromosomes other than the chr2 with the HOCS inserted. This code also provides quantification/counts of fragments of the 1L, 1L-2L linker, and 3R regions. These counts are used to generate the kinetics plots.

subnuc_zoom_footprint.R
This code generates a zoomed in dot plot (without a 2d cross correlation trace) of the region immediately around the HO cut site near PHO5. This is supplemental igure 6 of the paper.

genome_wide_genic_nucleosome_heatmap.R
This code generates a genome wide heatmap of aggregate nucleosome fuzziness within each gene body for every time point in the WT experiment. This is supplemental figure 7 of the paper.

nucleosome_pos_occupancy_plot.R
This code generates a plot that compares the occupancy of nucleosomes around the PHO5 dsb as well as the change in positioning of these nucleosomes between the pre-induction and 120min post induction samples. This is supplemental figure 8 of the paper.

facs_plots.R
This code will plot the .fcs files contained in raw_facs_files to generate supplemental figure 11 of the paper.
