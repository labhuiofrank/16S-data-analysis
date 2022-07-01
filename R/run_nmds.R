run_nmds <- function(mg, seed=seed, trymax=1000, dist="bray") {
  nmds_obj <- ordinate(mg, method="NMDS", distance=dist, trymax=trymax, seed=seed)
  components <- scores(nmds_obj)
  data <- sample_data(mg)[rownames(components),]
  data[, sprintf("NMDS%s", 1:2)] <- components
  stress <- nmds_obj$stress
  fi <- deparse(substitute(mg))
  file <- paste(fi,dist, sep="_")
  filename <- paste(file, "nmds.csv", sep="_")
  filepath <- paste(outdir,filename, sep="/")
  write.csv(data, file=filepath)
  print("NMDS stress=")
  print(stress)
  return(data)
}
