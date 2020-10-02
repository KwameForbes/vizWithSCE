integrateWithSingleCell <- function(res, dds) {
  # function written by Kwame Forbes, to assist with integration of
  # bulk DE results with Bioconductor single-cell RNA-seq datasets.
  # provides a menu of dataset options from pre-processed scRNA-seq
  # datasets. Downstream visualization is performed using the
  # vizWithSCE package:
  # https://github.com/KwameForbes/vizWithSCE
  stopifnot(is(res, "DESeqResults"))
  stopifnot(is(dds, "DESeqDataSet"))
  # figure out organism from 'dds' or 'res', either using tximeta metadata
  # or guesswork on the gene names in the results table
  tximetaOrg <- metadata(dds)$txomeInfo$organism
  if (!is.null(tximetaOrg)) {
    org <- if (tximetaOrg == "Homo sapiens") {
      "human"
    } else if (tximetaOrg == "Mus musculus") {
      "mouse"
    } else {
      stop("Only human and mouse are currently supported")
    }
  } else {
    test.gene <- rownames(res)[1]
    org <- if (substr(test.gene, 1, 4) == "ENSG") {
      "human"
    } else if (substr(test.gene, 1, 7) == "ENSMUSG") {
      "mouse"
    } else {
      stop("Only human and mouse are currently supported")
    }
  }
  message(paste("Your dataset appears to be", org))
  # read in table of options from vizWithSCE CSV file
  csv.file <- read.csv("R/singleCellTab.csv")
  #tab <- read.csv("csv.file")
  tab <- csv.file
  message(paste("Choose a",org,"single-cell dataset to integrate with (0 to cancel)"))
  tab <- tab[tab$org == org,]
  tab2 <- tab[,c("pkg","func","data", "pub","nCells")]
  tab2$data <- ifelse(is.na(tab2$data), "", tab2$data)
  rownames(tab2) <- seq_len(nrow(tab2))
  # print the table
  print(tab2)
  menuOpts <- ifelse(is.na(tab2$data), tab2$func, paste(tab2$func, tab2$data, sep="-"))
  ans <- menu(menuOpts)
  if (ans == 0) stop("No scRNA-seq dataset selected")
  pkg <- tab$pkg[ans]
  if (!requireNamespace(package=pkg, quietly=TRUE)) {
    message(paste0("Package: '",pkg, "' not installed."))
    ask <- askYesNo("Would you like to install package?")
    if (ask) {
      if (!requireNamespace(package="BiocManager", quietly=TRUE)) {
        stop("BiocManager required to install packages, install from CRAN")
      }
      BiocManager::install(pkg)
    } else {
      stop("Package would need to be installed")
    }
    if (!requireNamespace(package=pkg, quietly=TRUE)) {
      message("Package installed successfully")
    } else {
      stop("Package needs to be installed")
    }
  }
  # load package
  require(pkg, character.only=TRUE)
  # if the dataset is in the scRNAseq package...
  if (pkg == "scRNAseq") {
    # if only one dataset within the function...
    if (is.na(tab$data[ans])) {
      sce <- do.call(tab$func[ans], list(ensembl=TRUE))
    } else {
      sce <- do.call(tab$func[ans], list(which=tab$data[ans], ensembl=TRUE))
    }
  } else {
    # if only one dataset within the function...
    if (is.na(tab$data[ans])) {
      sce <- do.call(tab$func[ans], list())
    } else {
      sce <- do.call(tab$func[ans], list(dataset=tab$data[ans]))
    }
  }
  # return the original two objects and the SingleCellExperiment
  return(list(res=res, dds=dds, sce=sce))
}
getwd()

