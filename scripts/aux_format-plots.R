#==============================================================================
# Packages
#==============================================================================

load_packages <- function(package_list){
  for (p in package_list){
    suppressMessages(suppressWarnings(library(p, character.only = TRUE)))
  }
}

packages_default <- c("tidyverse", "cowplot", "patchwork", "RColorBrewer",
                      "ggpubr")

load_packages(packages_default)

#==============================================================================
# Auxiliary functions
#==============================================================================

adjust_lightness <- function(palette, lightnesses){
  # Adjust the lightnesses of a colour palette, while keeping hue and
  # saturation unchanged
  palette_hsl <- plotwidgets::col2hsl(palette)
  palette_hsl[3,] <- lightnesses
  return(plotwidgets::hsl2col(palette_hsl))
}

save_fig <- function(path, plot, width, aspect){
  cowplot::save_plot(path, plot, ncol = 1, nrow = 1, base_width = width,
                     base_height = width * aspect)
}

#==============================================================================
# Plotting theme information
#==============================================================================

# Import standard themes
source("scripts/aux_plot-theme.R")

#==============================================================================
# Palettes
#==============================================================================

#------------------------------------------------------------------------------
# Primary palette (Set2)
#------------------------------------------------------------------------------

# Colour names
palette_primary_names <- c("green1", "orange", "blue", "pink", "green2",
                           "yellow", "beige", "grey")

# Raw colour codes
palette_primary <- setNames(brewer.pal(8, "Set2"), palette_primary_names)
palette_primary_pale   <- setNames(brewer.pal(8, "Pastel2"),
                                   palette_primary_names)
palette_primary_dark   <- setNames(brewer.pal(8, "Dark2"),
                                   palette_primary_names)

# Fix pastel grey (too dark to read text)
palette_primary_pale["grey"] <- theme_base$legend.box.background$fill

# Intermediate blue & grey values
blue_mid <- adjust_lightness(palette_primary["blue"], 0.8)
grey_mid <- adjust_lightness(palette_primary["grey"], 0.8)
palette_primary_mid <- setNames(c(blue_mid, grey_mid),
                                c("blue", "grey"))

# Grey intermediate between normal and dark
grey_darkmid <- adjust_lightness(palette_primary["grey"], 0.6)
palette_primary_darkmid <- setNames(c(grey_darkmid), c("grey"))

# Darker normal blue (to distinguish from mid-blue)
palette_primary["blue"] <- adjust_lightness(palette_primary["blue"], 0.62)

#------------------------------------------------------------------------------
# Secondary palette (Set1)
#------------------------------------------------------------------------------

# Colour names
palette_secondary_names <- c("red", "blue", "green", "purple", "orange",
                           "yellow", "brown", "pink", "grey")

# Preferred lightnesses
palette_secondary_lightness <- c(0.63, 0.63, 0.63, 0.63, 
                                 0.73, 0.45, 0.63, 0.73, 0.63)
palette_secondary_darkness <- c(0.53, 0.53, 0.53, 0.53,
                                0.63, 0.35, 0.53, 0.63, 0.53)

# Raw colour codes
palette_secondary <- brewer.pal(9, "Set1") %>%
  adjust_lightness(palette_secondary_lightness) %>%
  setNames(palette_secondary_names)
palette_secondary_dark <- palette_secondary %>% 
  adjust_lightness(palette_secondary_darkness) %>%
  setNames(palette_secondary_names)