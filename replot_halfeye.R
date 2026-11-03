# Regenerate half-eye figures from saved plot_results.rds without rerunning simulations.
# Run from the directory containing the model scripts.
# Example: source("replot_halfeye.R"); replot_halfeye("results")
replot_halfeye <- function(root="results") {
  library(data.table)
  files <- list.files(root,pattern="^plot_results[.]rds$",recursive=TRUE,full.names=TRUE)
  if(!length(files)) stop("No plot_results.rds found under: ",root)
  e <- new.env(parent=globalenv())
  for(x in parse("simulation_model1.R"))
    if(is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("plot_halfeye"))) eval(x,e)
  for(file in files) {
    out <- dirname(file)
    model <- if(file.exists(file.path(out,"ancova_sim1_1.csv"))) 1L else if(file.exists(file.path(out,"ancova_2_1.csv"))) 2L else stop("Cannot identify model: ",out)
    f <- readRDS(file)$frequency
    for(sn in unique(f$scenario)) {
      g <- f[scenario==sn,.(n=sum(n)),by=.(method=sub("^bias_","",estimator),value)]
      name <- if(model==1L) switch(sn,general="new_sim1",common="new_sim2",ignorability="new_sim3") else
        switch(sn,general="new_sim2_1",common="new_sim2_2",ignorability="new_sim2_3")
      p <- e$plot_halfeye(g,name,out)
      stopifnot(identical(p$scales$get_scales("x")$limits,rev(c("Naive","ANCOVA","Change (F)","Change (A)","Change (D)"))))
    }
    message("Updated half-eye figures: ",out)
  }
}
