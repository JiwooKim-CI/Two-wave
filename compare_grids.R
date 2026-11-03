compare_grids <- function(root="results") {
  library(data.table)
  tables <- list()
  for(points in c(7L,9L)) for(model in 1:2) {
    path <- file.path(root,paste0("grid",points),paste0("model",model),"bias_summary.csv")
    if(!file.exists(path)) stop("Complete both grid analyses first: ",path)
    d <- fread(path)
    d[,`:=`(grid=points,model=paste0("model",model))]
    tables[[length(tables)+1L]]<-d
  }
  x<-rbindlist(tables)
  x[,rank_mean_abs_bias:=frank(mean_abs_bias,ties.method="min"),by=.(model,scenario,grid)]
  a<-x[grid==7,.(model,scenario,estimator,N7=N,MAB7=mean_abs_bias,rank7=rank_mean_abs_bias,max7=max_abs)]
  b<-x[grid==9,.(model,scenario,estimator,N9=N,MAB9=mean_abs_bias,rank9=rank_mean_abs_bias,max9=max_abs)]
  comparison<-merge(a,b,by=c("model","scenario","estimator"))
  comparison[,`:=`(MAB_difference=MAB9-MAB7,same_rank=rank7==rank9)]
  out<-file.path(root,"grid_comparison");dir.create(out,recursive=TRUE,showWarnings=FALSE)
  fwrite(comparison,file.path(out,"comparison_7_vs_9.csv"))
  checks<-x[(scenario=="common" & estimator=="bias_change_1")|(scenario=="ignorability" & estimator=="bias_ancova")]
  checks[,zero_bias_check:=max_abs<=1e-8]
  fwrite(checks,file.path(out,"identification_checks.csv"))
  p<-ggplot2::ggplot(x,ggplot2::aes(factor(grid),mean_abs_bias,colour=estimator,group=estimator))+
    ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::facet_grid(model~scenario,scales="free_y")+
    ggplot2::labs(x="Coefficient grid",y="Mean absolute population bias",colour="Estimator")+
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(out,"comparison_7_vs_9.pdf"),p,device=grDevices::pdf,width=11,height=6,useDingbats=FALSE)
  if(requireNamespace("ragg",quietly=TRUE))ggplot2::ggsave(file.path(out,"comparison_7_vs_9.png"),p,device=ragg::agg_png,width=11,height=6,dpi=300)
  print(comparison)
  invisible(comparison)
}
