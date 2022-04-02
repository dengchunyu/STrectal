#options(future.globals.maxSize = 10000 * 1024^2)
file_pre="/share/pub/dengcy/BRCA_organoid/RawData/Matrix/"
setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/Data_input_progress.r")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/Functions.r")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/library.r")
output_pre="/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress/"
setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")
load("BRCA_merged.RData")
plan("multiprocess", workers = N_WORKERS)
filename <- "./all_retransform_markers.rds"
if (!file.exists(filename)) {
  Allmarkers <- parallelFindAllMarkers(BRCA_merged)
  saveRDS(Allmarkers, filename)
} else {
  Allmarkers <- readRDS(filename)
}
