## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	collapse = TRUE,
	comment = "#>"
)

## ----include = FALSE-----------------------------------------------------
# install package
library(g3viz)

## ---- include = TRUE-----------------------------------------------------
# ============================================
# System file
maf.file <- system.file("extdata", "TCGA.BRCA.varscan.somatic.maf.gz", package = "g3viz")
# ============================================
# Read in MAF file
#   In addition to read data in, the built-in "readMAF" function also
#     1. parses "Mutation_Class" information according to the "Variant_Classification" column 
#        (also named "Mutation_Class" column in some files)
#     2. parses "AA_position" (amino-acid position) according to "HGVSp_Short" information
#        (e.g., p.V122Dfs*26, p.Q136P; also name "amino_acid_change" in some files)
mutation.dat <- readMAF(maf.file)

## ---- include = TRUE-----------------------------------------------------
# ============================================
# Chart 1
# use "default"" chart theme
# 
# Some features to check
#   1. zoom-in/out, pan, and tooltip tools in the lollipop panel
#   2. Use brush to select a reigon in the domain panel
#   3. click on lollipop pops to highlight findings
#   4. use legend to turn on/off certain mutation classes
#   5. click output buttons to save charts in SVG or PNG file
chart.options <- g3Lollipop.theme(theme.name = "default",
                                  title.text = "PIK3CA gene (default theme)")

g3Lollipop(mutation.dat,
           gene.symbol = "PIK3CA",
           plot.options = chart.options,
           output.filename = "default_theme")

## ---- include = TRUE-----------------------------------------------------
# ============================================
# Chart 2
# Visualize data according to "Variant_Classification" information
# use "nature" chart theme
chart.options <- g3Lollipop.theme(theme.name = "nature",
                                  title.text = "TP53 gene (nature theme)")

g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           factor.col = "Variant_Classification",
           # built-in gray-style output buttons
           btn.style = "gray",
           plot.options = chart.options,
           output.filename = "nature_theme")

## ---- include=TRUE-------------------------------------------------------
# ============================================
# load data
mutation.csv <- system.file("extdata", "ccle.csv", package = "g3viz")
# ============================================
# read in data
# specify "gene.symbol.col"    : column of gene symbol
#         "variant.class.col"  : column of variant class
#         "protein.change.col" : colum of protein change column
mutation.dat <- readMAF(mutation.csv,
                        gene.symbol.col = "Hugo_Symbol",
                        variant.class.col = "Variant_Classification",
                        protein.change.col = "amino_acid_change",
                        sep = ",")  # separator of csv file
# ============================================
# set up chart options
plot.options <- g3Lollipop.options(
  # chart options
  chart.width = 600,
  chart.type = "pie",
  chart.margin = list(left = 30, right = 20, top = 20, bottom = 30),
  chart.background = "#d3d3d3",
  # transition time
  transition.time = 300,
  # y-axis label
  y.axis.label = "# of TP53 gene mutations",
  axis.label.color = "#303030",
  axis.label.alignment = "end",
  axis.label.font = "italic 12px Serif",
  axis.label.dy = "-1.5em",
  # in-chart tick lines
  y.axis.line.color = "#303030",
  y.axis.line.width = 0.5,
  y.axis.line.style = "line",
  # legend
  legend = TRUE,
  legend.margin = list(left=20, right = 0, top = 10, bottom = 5),
  legend.interactive = TRUE,
  legend.title = "Variant classification",
  # lollipop track
  lollipop.track.height = 200,
  lollipop.track.background = "#d3d3d3",
  # lollipop pop
  lollipop.pop.min.size = 1,
  lollipop.pop.max.size = 8,
  # lollipop in-pop information
  lollipop.pop.info.limit = 5.5,
  lollipop.pop.info.dy = "0.24em",
  lollipop.pop.info.color = "white",
  # lollipop line
  lollipop.line.color = "#a9A9A9",
  lollipop.line.width = 3,
  # lollipop circle
  lollipop.circle.color = "#ffdead",
  lollipop.circle.width = 0.4,
  # lollipop click-on-pop highlight label
  lollipop.label.ratio = 2,
  lollipop.label.min.font.size = 12,
  # lollipop color scheme
  lollipop.color.scheme = "dark2",
  # highlight text angle
  highlight.text.angle = 60,
  # chart title
  title.color = "#303030",
  title.text = "TP53 gene (customized chart options)",
  title.font = "bold 12px monospace",
  title.alignment = "start",
  # annotation track
  anno.height = 16,
  anno.margin = list(top = 0, bottom = 0),
  anno.background = "#d3d3d3",
  # annotation track bar
  anno.bar.fill = "#a9a9a9",
  anno.bar.margin = list(top = 4, bottom = 4),
  # protein domain
  domain.color.scheme = "pie5",
  domain.margin = list(top = 2, bottom = 2),
  domain.text.color = "white",
  domain.text.font = "italic 8px Serif",
  # selection brush
  brush = TRUE,
  brush.selection.background = "#F8F8FF",
  brush.selection.opacity = 0.3,
  brush.border.color = "#a9a9a9",
  brush.border.width = 1,
  brush.handler.color = "#303030",
  # tooltip
  tooltip = TRUE,
  # zoom
  zoom = TRUE
)

g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           protein.change.col = "amino_acid_change",
           # built-in blue-style output buttons
           btn.style = "blue",
           plot.options = plot.options,
           output.filename = "customized_plot")

## ---- include=TRUE-------------------------------------------------------
# ============================================
# get genomic mutation data of "msk_impact_2017"" study from cBioPortal
mutation.dat <- getMutationsFromCbioportal("msk_impact_2017", "TP53")

# g3viz has built-in "cbioportal" theme
plot.options <- g3Lollipop.theme(theme.name = "cbioportal",
                                 title.text = "TP53 gene (cbioportal theme)",
                                 y.axis.label = "# of TP53 Mutations")

g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           btn.style = "blue",
           plot.options = plot.options,
           output.filename = "cbioportal_theme")

## ---- include = TRUE-----------------------------------------------------
# ===============================
# TP54 has single UniProt entry
hgnc2pfam("TP53", output.format = "list")
# ===============================
# GNAS gene has multiple UniProt entries
#   `guess = TRUE`: the Pfam domain compositon of the longest 
#                   UniProt protein is returned
hgnc2pfam("GNAS", guess = TRUE)

## ---- out.width = "620px", echo = FALSE----------------------------------
knitr::include_graphics("figures/color_scheme.png")

## ---- include = TRUE-----------------------------------------------------
# ============================================
# Load system sample data
maf.file <- system.file("extdata", "TCGA.BRCA.varscan.somatic.maf.gz", package = "g3viz")
mutation.dat <- readMAF(maf.file)

## ---- include = TRUE-----------------------------------------------------
# ============================================
# simple theme
chart.options <- g3Lollipop.theme(theme.name = "simple",
                                  title.text = "TP53 gene (simple theme)")

g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           plot.options = chart.options,
           btn.style = "gray",
           output.filename = "simple_theme")

## ---- include = TRUE-----------------------------------------------------
# ============================================
# ggplot2 theme
chart.options <- g3Lollipop.theme(theme.name = "ggplot2",
                                  title.text = "TP53 gene (ggplot2 theme)")

g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           plot.options = chart.options,
           btn.style = "gray",
           output.filename = "ggplot2_theme")

## ---- include=TRUE-------------------------------------------------------
sessionInfo()

