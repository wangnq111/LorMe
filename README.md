<!-- badges: start -->
[![CRAN
Version](https://www.r-pkg.org/badges/version/LorMe)](https://cran.r-project.org/package=LorMe)
<!-- badges: end -->
# LorMe
LorMe package for R：Lightening One-code Resolving Microbial Ecology Solution. LorMe keeps updating and maintanence. 

# Suggestions and bug report
2434066068@qq.com

# Installation
Standard version from CRAN
```{R}
install.packages("LorMe")
```
Development version from GitHub :stuck_out_tongue_closed_eyes:
```{R}
if (!require(remotes)) install.packages("remotes")
remotes::install_github("wangnq111/LorMe")
```


# Illustration
[Getting Started](https://wangnq111.github.io/Gettingstarted.html)
## Chinese version
Chinese illustration available at 🏮 [LorMe中文版说明书](https://rural-dianella-be0.notion.site/LorMe-aac2ba66a3bf46bd89c103e78550e6f4) 🏮
## Visualization

### Community Feature
#### Alpha diversity
<img src="./man/figures/README_Alpha.png" width="60%" style="display: block; margin: auto;" />

#### Community Structure
Support three styles: ellipse (as in the below),stick and polygon
<img src="./man/figures/README_structure.png" width="60%" style="display: block; margin: auto;" />

#### Community Composition
Support three styles: Bar plot (as in the below), Area plot and alluvial plot
<img src="./man/figures/README_community.png" width="60%" style="display: block; margin: auto;" />

### Differential Analysis
#### Differential Bar
<img src="./man/figures/README_Diff.png" width="70%" style="display: block; margin: auto;" />

#### Volcano Plot
Support Fold change-FDR plot and Mean-Fold change plot
<img src="./man/figures/README_volcano.png" width="60%" style="display: block; margin: auto;" />

#### Manhatton Plot
Support both classical style and circular style
<img src="./man/figures/README_Manhattan.png" width="100%" style="display: block; margin: auto;" />

### Network analysis
#### Classical network
Painted top five largest modules
<img src="./man/figures/README_network.png" width="60%" style="display: block; margin: auto;" />

#### Meta network
Painted differential taxon
<img src="./man/figures/README_Metanetwork.png" width="60%" style="display: block; margin: auto;" />

#### Module composition pie chart
<img src="./man/figures/README_Modulepie.png" width="60%" style="display: block; margin: auto;" />
