<img src="./krigzy.png" align = "right" height=250/>

# Krigzy

A friendly Shiny application to learn ordinary kriging.

To use Krigzy you need to have following packages installed: <a href="https://github.com/rstudio/shiny">shiny</a>, <a href="https://github.com/dreamRs/shinyWidgets">shinyWidgets</a>, <a href="https://github.com/daattali/shinyalert/">shinyalert</a>, <a href="https://github.com/r-spatial/gstat">gstat</a>, <a href="https://github.com/r-spatial/stars/">stars</a>, <a href="https://github.com/r-spatial/sf">sf</a>, <a href="https://github.com/jrowen/rhandsontable">rhandsontable</a>, <a href="https://github.com/tidyverse/ggplot2">ggplot2</a>, <a href="https://github.com/cran/plotrix">plotrix</a>.

### Prerequsites:

```r
install.packages(c("shiny", "shinyWidgets",
                   "shinyalert", "gstat",
                   "stars", "sf", "rhandsontable", 
                   "ggplot2", "plotrix"))

```

### Download and run:

```r
library(shiny)
runGitHub("Krigzy", "themercerus", subdir = "Krigzy/", ref = "main")
```
