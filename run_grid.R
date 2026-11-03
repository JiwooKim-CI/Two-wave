# Run from the directory containing the two model scripts.
run_grid <- function(points, root="results", draw_plots=TRUE) {
  stopifnot(points %in% c(3L,7L,9L))
  old <- Sys.getenv(c("N_POINTS","OUTPUT_DIR","DRAW_PLOTS"),unset=NA_character_)
  on.exit(for(n in names(old)) if(is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv,setNames(list(old[[n]]),n)))
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  for(model in 1:2) {
    out <- file.path(root,paste0("grid",points),paste0("model",model))
    Sys.setenv(N_POINTS=as.character(points),OUTPUT_DIR=out,DRAW_PLOTS=if(draw_plots)"TRUE" else "FALSE")
    message("Running model ",model," on the ",points,"-point grid: ",out)
    status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(paste0("simulation_model",model,".R"))))
    if(status!=0L) stop("Analysis failed; rerun this command to resume completed beta blocks.")
  }
}
