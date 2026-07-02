# Lake-Chamo-Network-Resilience
# Lake Chamo Network Resilience Framework

R code for the paper: 
*"Network Destabilization and Ecological Memory as Early Warning Signals for Algal Blooms in a Hypereutrophic Tropical Lake: A Framework for Lake Chamo, Ethiopia"*

---

## Authors
- Habtamu Tadesse
- Kassahun Atalaya
- Gemamaw Tesfaye
- Gashaw Tesfaye

**Affiliation:** Arba Minch University, Ethiopia

---

## Repository Contents

| File | Description |
|------|-------------|
| `lake_chamo_analysis.R` | Complete R code for network analysis, keystone species identification, eco-network energy calculation, and satellite NDVI validation |

---

## Requirements

- **R version:** 4.0 or higher
- **Required R packages:**
  - `igraph` - Network construction and analysis
  - `terra` - Raster processing
  - `sf` - Spatial data handling
  - `ggplot2` - Visualization
  - `tidyterra` - ggplot2 for rasters
  - `rstac` - Access satellite data
  - `dplyr` - Data manipulation
  - `tidyr` - Data reshaping
  - `patchwork` - Combining plots
  - `RColorBrewer` - Color palettes
  - `ggspatial` - Scale bars and north arrows
  - `Hmisc` - Correlation matrices

To install all packages at once:
```r
install.packages(c("igraph", "terra", "sf", "ggplot2", "tidyterra", 
                   "rstac", "dplyr", "tidyr", "patchwork", 
                   "RColorBrewer", "ggspatial", "Hmisc"))
