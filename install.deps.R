install.packages(
  pkgs=c(
    "http://cran.us.r-project.org/src/contrib/Matrix_1.1-4.tar.gz",
    "http://cran.us.r-project.org/src/contrib/Rcpp_0.11.2.tar.gz"
  ),
  repos=NULL,
  lib="~/R_libs"
)
install.packages(
  c(
    "plyr", "caret","corrplot",
    "data.table"
  ),
  repos="http://cran.us.r-project.org/",
  lib="~/R_libs"
)
