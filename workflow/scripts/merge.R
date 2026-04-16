## collect QDNAseq objects or data frames for single cell profiles into single object
require(BiocGenerics, quietly=TRUE)
require(gtools, quietly=TRUE)
require(dplyr, quietly=TRUE)
if (interactive()){
  BASEDIR="~/scAbsolute"
} else {
  BASEDIR="/opt/scAbsolute"
}
source(file.path(BASEDIR, "R/core.R"))
args = commandArgs(trailingOnly=TRUE)
indexSort = TRUE
options(future.globals.maxSize= 4096*1024^2)
rdsFiles = tail(args, n=-1)
if (length(args) >= 3){
  d = args[2]
  f = args[3]  
}else{
  d = args[2]
  f = NULL
}
if (dir.exists(d) && (is.null(f) || !file.exists(f))){
  setwd(d)
  rdsFiles = list.files(pattern=paste0(f, ".rds$"), recursive=FALSE)
}else{
  stopifnot(all(endsWith(rdsFiles, ".rds")))
  stopifnot(all(file.exists(rdsFiles)))
  if(indexSort){
    rdsFiles = sort(rdsFiles, decreasing=FALSE)
  }
  filesToRemove <- character(0)
  for (file in rdsFiles) {
    if (file.info(file)$size == 0) {
      file.remove(file)
      filesToRemove <- c(filesToRemove, file)
      cat("File removed:", file, "\n")
    }
  }
  rdsFiles <- setdiff(rdsFiles, filesToRemove)
}
if (is.data.frame(readRDS(rdsFiles[[1]]))) {
  rdsData = base::lapply(rdsFiles, readRDS)
  mergedRDS = do.call(rbind, rdsData)
} else {
  rdsData = base::sapply(rdsFiles, readRDS, USE.NAMES = FALSE)
  rdsData = unlist(rdsData)

  # Normalise pData schemas across all cells — union all columns, fill missing with NA
  all_cols <- unique(unlist(lapply(rdsData, function(x) colnames(Biobase::pData(x)))))
  ref_cols <- colnames(Biobase::pData(rdsData[[1]]))
  # Use column order from a passing cell (fit_flag=TRUE) as reference if possible
  passing <- Filter(function(x) isTRUE(Biobase::pData(x)[["fit_flag"]]), rdsData)
  if (length(passing) > 0) ref_cols <- colnames(Biobase::pData(passing[[1]]))
  extra_cols <- setdiff(all_cols, ref_cols)
  final_cols <- c(ref_cols, extra_cols)

  rdsData <- lapply(rdsData, function(x) {
    desc <- Biobase::pData(x)
    missing <- setdiff(final_cols, colnames(desc))
    for (col in missing) desc[[col]] <- NA
    Biobase::pData(x) <- desc[, final_cols, drop=FALSE]
    x
  })


# Filter out failed cells (fit_flag=FALSE) BEFORE combining
# so combineQDNASets receives a consistent set of passing cells
failed_cells_df <- data.frame(name=character(), failure_reason=character(), rpc=numeric())
if (!is.data.frame(rdsData[[1]])) {
  fit_flags <- sapply(rdsData, function(x) isTRUE(Biobase::pData(x)[["fit_flag"]]))
  if (any(!fit_flags)) {
    failed <- rdsData[!fit_flags]
    failed_cells_df <- do.call(rbind, lapply(failed, function(x) {
      pd <- Biobase::pData(x)
      data.frame(name=pd[["name"]], failure_reason=as.character(pd[["failure_reason"]]),
                 rpc=pd[["rpc"]], stringsAsFactors=FALSE)
    }))
    cat("Excluding", nrow(failed_cells_df), "failed cells (fit_flag=FALSE) from merged RDS\n")
    rdsData <- rdsData[which(fit_flags)]
  }
}
filesToRemove <- c(filesToRemove, failed_cells_df$name)

mergedRDS = combineQDNASets(rdsData)
}
write.csv(if(nrow(failed_cells_df)>0) failed_cells_df else filesToRemove,
  file=paste0(dirname(args[1]),"/failed_cells.csv"), row.names=FALSE)
# Drop columns that are all-NA across all remaining cells (artefacts from schema normalisation)
if (!is.data.frame(mergedRDS)) {
  pd <- Biobase::pData(mergedRDS)
  all_na_cols <- names(which(sapply(pd, function(col) all(is.na(col)))))
  if (length(all_na_cols) > 0) {
    cat("Dropping", length(all_na_cols), "all-NA columns from merged RDS\n")
    Biobase::pData(mergedRDS) <- pd[, !colnames(pd) %in% all_na_cols, drop=FALSE]
  }
}
saveRDS(mergedRDS, file=args[1])
