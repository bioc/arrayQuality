.First.lib <- function(libname, pkgname) {

  print(pkgname)
  print("Checking for missing packages")

  packlist <- installed.packages()
  
  ## Check if required BioConductor packages are installed

  requiredPacks <- c("marray", "limma", "hexbin")
  missingPacks <- requiredPacks[!requiredPacks %in% packlist]

  if(length(missingPacks) > 0)
    {
      print("Downloading missing packages from BioConductor")
      if(! "reposTools" %in% packlist)
        {
          #Use getBioC
          source("http://www.bioconductor.org/getBioC.R")
          getBioC(libName="reposTools", bundle=FALSE, getAllDeps=TRUE)
          install.packages2(pkgs=missingPacks, getAllDeps=TRUE)
        }  
      else{
        require(reposTools)
        install.packages2(pkgs=missingPacks, getAllDeps=TRUE)
      }
    }
  
  #Test for packages from CRAN
  if(! "mclust" %in% packlist)
    {
      print("Downloading missing packages from CRAN")
      install.packages("mclust")
    }
  
  rm(packlist, requiredPacks, missingPacks); gc(TRUE)
  print("Loading required packages")
  require("marray", quietly=TRUE)
  require("limma", quietly=TRUE)
  require("mclust", quietly=TRUE)
  require("hexbin", quietly=TRUE)
}


