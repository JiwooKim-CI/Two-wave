
# Install once: install.packages(c("data.table", "ggplot2", "ggdist", "jtools", "ragg"))
library(data.table)
N_POINTS <- as.integer(Sys.getenv("N_POINTS", "7"))
INNER_DIM <- as.integer(Sys.getenv("INNER_DIM", "6"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", paste0("model1_corrected_", format(Sys.time(), "%Y%m%d_%H%M%S")))
DRAW_PLOTS <- Sys.getenv("DRAW_PLOTS", "TRUE") == "TRUE"
setDTthreads(as.integer(Sys.getenv("R_THREADS", "1")))
stopifnot(N_POINTS %in% c(3L, 7L, 9L), INNER_DIM >= 1, INNER_DIM <= 6)
CHECKPOINT <- file.path(OUTPUT_DIR,"checkpoint.rds")
if (dir.exists(OUTPUT_DIR) && !file.exists(CHECKPOINT)) stop("Existing directory has no checkpoint; choose a new OUTPUT_DIR.")
if (DRAW_PLOTS) for (pkg in c("ggplot2","ggdist","jtools"))
  if (!requireNamespace(pkg, quietly=TRUE)) stop("Install package: ",pkg)
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)
# Three points are for testing only.
ng <- if (N_POINTS==3L) c(-.1,0,.1) else {
  positive <- .1^(1/2^(0:((N_POINTS-3L)/2L)))
  c(-positive,0,rev(positive))
}
p_m2_x <- .2
p_y2_x <- .2
bias_names <- c("bias_ancova","bias_naive","bias_change","bias_change_1","bias_change_2")


calculate_chunk <- function(grid1) {
grid1[, var_eps_m1 := 1 - p_m1_u^2]
grid1[, var_eps_y1 := 1 - p_y1_u^2 - p_y1_m1^2 - 2 * p_y1_u * p_m1_u * p_y1_m1]
grid1[, var_eps_m2 := (1 - p_m2_m1^2 - p_m2_u^2 - p_m2_x^2 - p_m2_y1^2 + (1 - p_m2_m1^2 - p_m2_u^2 - p_m2_x^2 - p_m2_y1^2 - 
    p_m2_m1 * (2 * p_m2_u * p_m1_u + 2 * p_m2_y1 * p_y1_m1 + 2 * p_m2_y1 * p_y1_u * p_m1_u) - p_m2_u * 
    (2 * p_m2_m1 * p_m1_u + 2 * p_m2_y1 * p_y1_m1 * p_m1_u + 2 * p_m2_y1 * p_y1_u) - p_m2_y1 * (2 * p_m2_m1 * 
    p_y1_m1 + 2 * p_m2_m1 * p_m1_u * p_y1_u + 2 * p_m2_u * p_m1_u * p_y1_m1 + 2 * p_m2_u * p_y1_u)))/2]
grid1[, var_eps_y2 := (1 - p_y2_y1^2 - p_y2_m2^2 - p_y2_u^2 - p_y2_x^2 - p_y2_m1^2 + (1 - p_y2_y1^2 - p_y2_m2^2 - p_y2_u^2 - 
    p_y2_x^2 - p_y2_m1^2 - p_y2_m2 * (2 * p_y2_m1 * p_m2_m1 + 2 * p_y2_m1 * p_y1_m1 * p_m2_y1 + 2 * p_y2_m1 * 
    p_m1_u * p_m2_u + 2 * p_y2_m1 * p_m1_u * p_y1_u * p_m2_y1 + 2 * p_y2_u * p_m1_u * p_m2_m1 + 2 * p_y2_u * 
    p_m1_u * p_y1_m1 * p_m2_y1 + 2 * p_y2_u * p_m2_u + 2 * p_y2_u * p_y1_u * p_m2_y1 + 2 * p_y2_x * p_m2_x + 
    2 * p_y2_y1 * p_m2_y1 + 2 * p_y2_y1 * p_y1_m1 * p_m2_m1 + 2 * p_y2_y1 * p_y1_m1 * p_m1_u * p_m2_u + 
    2 * p_y2_y1 * p_y1_u * p_m1_u * p_m2_m1 + 2 * p_y2_y1 * p_y1_u * p_m2_u) - p_y2_m1 * (2 * p_y2_m2 * 
    p_m2_m1 + 2 * p_y2_m2 * p_m2_u * p_m1_u + 2 * p_y2_m2 * p_m2_y1 * p_y1_m1 + 2 * p_y2_m2 * p_m2_y1 * 
    p_y1_u * p_m1_u + 2 * p_y2_u * p_m1_u + 2 * p_y2_y1 * p_y1_m1 + 2 * p_y2_y1 * p_y1_u * p_m1_u) - 
    p_y2_y1 * (2 * p_y2_m1 * p_y1_m1 + 2 * p_y2_m1 * p_m1_u * p_y1_u + 2 * p_y2_m2 * p_m2_m1 * p_y1_m1 + 
        2 * p_y2_m2 * p_m2_m1 * p_m1_u * p_y1_u + 2 * p_y2_m2 * p_m2_u * p_m1_u * p_y1_m1 + 2 * p_y2_m2 * 
        p_m2_u * p_y1_u + 2 * p_y2_m2 * p_m2_y1 + 2 * p_y2_u * p_m1_u * p_y1_m1 + 2 * p_y2_u * p_y1_u) - 
    p_y2_u * (2 * p_y2_m1 * p_m1_u + 2 * p_y2_m2 * p_m2_m1 * p_m1_u + 2 * p_y2_m2 * p_m2_u + 2 * p_y2_m2 * 
        p_m2_y1 * p_y1_m1 * p_m1_u + 2 * p_y2_m2 * p_m2_y1 * p_y1_u + 2 * p_y2_y1 * p_y1_m1 * p_m1_u + 
        2 * p_y2_y1 * p_y1_u) - 2 * p_y2_m2 * p_m2_x * p_y2_x))/2]


grid1 <- grid1[var_eps_m1 > 1e-10 & var_eps_y1 > 1e-10 & var_eps_m2 > 1e-10 & var_eps_y2 > 1e-10]
for (v in c("m1","y1","m2","y2")) grid1[, (paste0("eps_",v)) := sqrt(get(paste0("var_eps_",v)))]

grid1[, rm2y2 := p_y2_m2 + p_m2_m1 * p_y1_m1 * p_y2_y1 + p_m2_m1 * p_y2_m1 + p_m2_m1 * p_m1_u * p_y1_u * p_y2_y1 + p_m2_m1 * 
    p_m1_u * p_y2_u + p_m2_u * p_m1_u * p_y1_m1 * p_y2_y1 + p_m2_u * p_m1_u * p_y2_m1 + p_m2_u * p_y1_u * 
    p_y2_y1 + p_m2_u * p_y2_u + p_m2_x * p_y2_x + p_m2_y1 * p_y2_y1 + p_m2_y1 * p_y1_m1 * p_y2_m1 + p_m2_y1 * 
    p_y1_m1 * p_m1_u * p_y2_u + p_m2_y1 * p_y1_u * p_m1_u * p_y2_m1 + p_m2_y1 * p_y1_u * p_y2_u]
grid1[, rm1m2 := p_m2_m1 + p_m1_u * p_m2_u + p_y1_m1 * p_m2_y1 + p_m1_u * p_y1_u * p_m2_y1]
grid1[, rm1y2 := p_m2_m1 * p_y2_m2 + p_y1_m1 * p_m2_y1 * p_y2_m2 + p_y1_m1 * p_y2_y1 + p_y2_m1 + p_m1_u * p_m2_u * p_y2_m2 + 
    p_m1_u * p_y1_u * p_m2_y1 * p_y2_m2 + p_m1_u * p_y1_u * p_y2_y1 + p_m1_u * p_y2_u]
grid1[, ry1m2 := p_m2_y1 + p_m2_m1 * p_y1_m1 + p_m2_m1 * p_m1_u * p_y1_u + p_m2_u * p_m1_u * p_y1_m1 + p_m2_u * p_y1_u]
grid1[, ry1y2 := p_y2_m2 * p_m2_y1 + p_y2_y1 + p_y2_m2 * p_m2_m1 * p_y1_m1 + p_y2_m1 * p_y1_m1 + p_y2_m2 * p_m2_m1 * p_m1_u * 
    p_y1_u + p_y2_m1 * p_m1_u * p_y1_u + p_y2_m2 * p_m2_u * p_m1_u * p_y1_m1 + p_y2_u * p_m1_u * p_y1_m1 + 
    p_y2_m2 * p_m2_u * p_y1_u + p_y2_u * p_y1_u]
grid1[, rm1y1 := p_y1_m1 + p_m1_u * p_y1_u]
grid1[, rm1x := 0]
grid1[, ry1x := 0]
grid1[, rm2x := p_m2_x]
grid1[, ry2x := p_m2_x * p_y2_m2 + p_y2_x]
grid1[, ry1m2_m1 := (ry1m2 - rm1y1 * rm1m2)/sqrt((1 - rm1y1^2) * (1 - rm1m2^2))]
grid1[, ry1y2_m1 := (ry1y2 - rm1y1 * rm1y2)/sqrt((1 - rm1y1^2) * (1 - rm1y2^2))]
grid1[, rm2y2_m1 := (rm2y2 - rm1m2 * rm1y2)/sqrt((1 - rm1m2^2) * (1 - rm1y2^2))]
grid1[, ry1x_m1 := (ry1x - rm1y1 * rm1x)/sqrt((1 - rm1y1^2) * (1 - rm1x^2))]
grid1[, rxy2_m1 := (ry2x - rm1x * rm1y2)/sqrt((1 - rm1x^2) * (1 - rm1y2^2))]
grid1[, rxm2_m1 := (rm2x - rm1x * rm1m2)/sqrt((1 - rm1x^2) * (1 - rm1m2^2))]
grid1[, sd_m2_m1y1 := sqrt(1 - (ry1m2^2 + rm1m2^2 - 2 * ry1m2 * rm1m2 * rm1y1)/(1 - rm1y1^2))]
grid1[, sd_y2_m1y1 := sqrt(1 - (ry1y2^2 + rm1y2^2 - 2 * ry1y2 * rm1y2 * rm1y1)/(1 - rm1y1^2))]
grid1[, by2m2_x := (rm2y2 - ry2x * rm2x)/((1 - rm2x^2))]
grid1[, rm2y2_m1y1 := (rm2y2_m1 - ry1m2_m1 * ry1y2_m1)/sqrt((1 - ry1m2_m1^2) * (1 - ry1y2_m1^2))]
grid1[, rxy2_m1y1 := (rxy2_m1 - ry1x_m1 * ry1y2_m1)/sqrt((1 - ry1x_m1^2) * (1 - ry1y2_m1^2))]
grid1[, rxm2_m1y1 := (rxm2_m1 - ry1x_m1 * ry1m2_m1)/sqrt((1 - ry1x_m1^2) * (1 - ry1m2_m1^2))]
grid1[, bm2y2_xm1y1 := ((rm2y2_m1y1 - rxm2_m1y1 * rxy2_m1y1)/(1 - rxm2_m1y1^2)) * sd_y2_m1y1/sd_m2_m1y1]


grid1 <- grid1[sd_m2_m1y1^2 * (1-rxm2_m1y1^2) > 1e-10]

# Variance of raw mediator change conditional on X
grid1[, var_dm_x :=
        2 - 2 * rm1m2 - (rm2x - rm1x)^2]

# Covariance between the two raw change scores conditional on X
grid1[, cov_dy_dm_x :=
        (rm2y2 - rm1y2 - ry1m2 + rm1y1) -
        (ry2x - ry1x) * (rm2x - rm1x)]

# Variance of M2 conditional on X
grid1[, var_m2_x :=
        1 - rm2x^2]

# Covariance between Delta Y and M2 conditional on X
grid1[, cov_dy_m2_x :=
        (rm2y2 - ry1m2) -
        (ry2x - ry1x) * rm2x]

# Remove parameter combinations producing invalid conditional variances
grid1 <- grid1[
  var_dm_x > 1e-10 &
  var_m2_x > 1e-10
]

grid1[,`:=`(change = (rm2x - rm1x) *
        cov_dy_dm_x / var_dm_x,
            
            change_1 = rm2x *
        cov_dy_m2_x / var_m2_x,
            ancova = p_m2_x * bm2y2_xm1y1 ,
            truth = p_m2_x * p_y2_m2,
            naive =p_m2_x * by2m2_x)]



# Convert to a tibble for better display if desired

grid1[, change_2 := 2 * change]

# computing bias
grid1[,`:=`(bias_change = change - truth, ### change(F)
            bias_change_2 = change_2 - truth, ### change(A)
            bias_ancova = ancova - truth,
            bias_naive = naive - truth,
            bias_change_1 = change_1 - truth)] ### change(D)
# add the column showing the degree of violation of the ignorability assumption
grid1$diff <- abs(grid1$p_m2_u * grid1$p_y2_u)


  grid1 <- grid1[is.finite(ancova) & is.finite(naive) & is.finite(change) & is.finite(change_1) & is.finite(change_2)]
  grid1
}

# Accumulate frequencies rather than retain billions of raw rows.
frequency <- summary_table <- wins <- curves <- NULL
counts <- list()
add_results <- function(grid1small1, scenario) {
  if (!nrow(grid1small1)) return(invisible(NULL))
  if (scenario=="common" && max(abs(grid1small1$bias_change_1))>1e-8) stop("Common-trend check failed.")
  if (scenario=="ignorability" && max(abs(grid1small1$bias_ancova))>1e-8) stop("Ignorability check failed.")
  for (v in bias_names) {
    b <- grid1small1[[v]]
    graph1 <- data.table(scenario=scenario, estimator=v, value=round(b,4), diff=round(grid1small1$diff,4))
    graph1 <- graph1[,.(n=as.double(.N)),by=.(scenario,estimator,value,diff)]
    frequency <<- rbindlist(list(frequency,graph1),use.names=TRUE)[,.(n=sum(n)),by=.(scenario,estimator,value,diff)]
    z <- data.table(scenario=scenario,estimator=v,N=as.double(length(b)),sum_bias=sum(b),sum_abs=sum(abs(b)),max_abs=max(abs(b)))
    summary_table <<- rbindlist(list(summary_table,z),use.names=TRUE)[,.(N=sum(N),sum_bias=sum(sum_bias),sum_abs=sum(sum_abs),max_abs=max(max_abs)),by=.(scenario,estimator)]
    if (scenario=="general") {
      z <- data.table(estimator=v,diff=round(grid1small1$diff,12),ab=abs(b))[,.(N=as.double(.N),sum_abs=sum(ab)),by=.(estimator,diff)]
      curves <<- rbindlist(list(curves,z),use.names=TRUE)[,.(N=sum(N),sum_abs=sum(sum_abs)),by=.(estimator,diff)]
    }
  }
  ab <- abs(as.matrix(grid1small1[,..bias_names]))
  tied <- abs(ab-apply(ab,1,min))<=1e-10
  grid1small1[,Min:=ifelse(rowSums(tied)>1,"Tie",bias_names[max.col(tied,ties.method="first")])]
  z <- grid1small1[,.(n=as.double(.N)),by=Min][,scenario:=scenario]
  wins <<- rbindlist(list(wins,z),use.names=TRUE)[,.(n=sum(n)),by=.(scenario,Min)]
}

# Resume completed beta blocks with the same model/grid settings.
config <- list(version="complete_v1", model=1L, points=N_POINTS, inner=INNER_DIM,
               calculation=deparse(body(calculate_chunk)), accumulation=deparse(body(add_results)))
resume <- NULL
if (file.exists(CHECKPOINT)) {
  resume <- readRDS(CHECKPOINT)
  if (!identical(resume$config,config)) stop("Checkpoint settings/code differ; choose a new OUTPUT_DIR.")
  frequency <- resume$frequency; summary_table <- resume$summary_table
  wins <- resume$wins; curves <- resume$curves; counts <- resume$counts
  message("Resuming after ",resume$scenario," beta index ",resume$beta)
}
save_checkpoint <- function(scenario,i,candidates,retained) {
  tmp <- paste0(CHECKPOINT,".tmp")
  saveRDS(list(config=config,scenario=scenario,beta=i,candidates=candidates,retained=retained,
    frequency=frequency,summary_table=summary_table,wins=wins,curves=curves,counts=counts),tmp)
  if (!file.rename(tmp,CHECKPOINT)) stop("Could not save checkpoint.")
}

# Preserve CJ enumeration; only split it into smaller grids.
for (scenario in c("general","common")) {
  free_names <- c("p_m2_m1","p_y2_y1","p_y2_u","p_m2_u","p_m2_y1","p_y1_u","p_y2_m1","p_m1_u","p_y1_m1")
  if (scenario=="common") free_names <- setdiff(free_names,c("p_y2_u","p_y2_m1","p_m2_y1"))
  inner_names <- tail(free_names,min(INNER_DIM,length(free_names)))
  outer_names <- setdiff(free_names,inner_names)
  outer_grid <- if (length(outer_names)) do.call(CJ,setNames(rep(list(ng),length(outer_names)),outer_names)) else data.table(dummy=1)
  if (!is.null(resume) && match(scenario,c("general","common")) < match(resume$scenario,c("general","common"))) next
  same_resume <- !is.null(resume) && scenario==resume$scenario
  candidates <- if(same_resume) resume$candidates else 0
  retained <- if(same_resume) resume$retained else 0
  for (i in seq_along(ng)) {
    if (same_resume && i<=resume$beta) next
    p_y2_m2 <- ng[i]
    for (j in seq_len(nrow(outer_grid))) {
      grid1 <- do.call(CJ,setNames(rep(list(ng),length(inner_names)),inner_names))
      for (v in outer_names) grid1[, (v):=outer_grid[[v]][j]]
      grid1[, `:=`(p_y2_m2=p_y2_m2,p_m2_x=p_m2_x,p_y2_x=p_y2_x)]
      if (scenario=="common") grid1[, `:=`(p_m2_y1=0,p_y2_u=(1-p_y2_y1)*p_y1_u,p_y2_m1=(1-p_y2_y1)*p_y1_m1)]
      candidates <- candidates+nrow(grid1)
      grid1 <- calculate_chunk(grid1)
      retained <- retained+nrow(grid1)
      add_results(grid1,scenario)
      if (scenario=="general") {
        grid1small1 <- grid1[p_m2_u*p_y2_u==0]
        add_results(grid1small1,"ignorability")
      }
    }
    counts[[scenario]] <- data.table(scenario=scenario,candidates=candidates,retained=retained)
    save_checkpoint(scenario,i,candidates,retained)
    message(scenario,": beta ",i,"/",length(ng)," completed; candidates=",candidates)
  }
  counts[[scenario]] <- data.table(scenario=scenario,candidates=candidates,retained=retained)
}
summary_table[, `:=`(mean_bias=sum_bias/N,mean_abs_bias=sum_abs/N)]
wins[,proportion:=n/sum(n),by=scenario]
curves[,mean_abs_bias:=sum_abs/N]
fwrite(summary_table,file.path(OUTPUT_DIR,"bias_summary.csv"))
fwrite(wins,file.path(OUTPUT_DIR,"best_estimator.csv"))
fwrite(rbindlist(counts),file.path(OUTPUT_DIR,"counts.csv"))
fwrite(curves,file.path(OUTPUT_DIR,"violation_curve.csv"))
saveRDS(list(frequency=frequency,summary=summary_table,curves=curves),file.path(OUTPUT_DIR,"plot_results.rds"))
# Keep the original per-estimator CSV naming scheme.
stems <- c(bias_ancova="ancova",bias_naive="naive",bias_change="change",bias_change_1="change1",bias_change_2="change2")
for (sn in unique(frequency$scenario)) for (v in bias_names) {
  graph1 <- frequency[scenario==sn & estimator==v,.(value,diff,n)]
  graph1[,perc:=n/sum(n)]
  suffix <- switch(sn,general="sim1_1",common="sim1_c",ignorability="sim1_i")
  fwrite(graph1,file.path(OUTPUT_DIR,paste0(stems[[v]],"_",suffix,".csv")))
}
plot_halfeye <- function(d, name, output_dir, save_tex = FALSE) {
  labels <- c(naive = "Naive", change = "Change (F)", change_2 = "Change (A)",
              change_1 = "Change (D)", ancova = "ANCOVA")
  display_order <- c("Naive", "ANCOVA", "Change (F)", "Change (A)", "Change (D)")
  # coord_flip displays the last factor level at the top.
  graph <- copy(d)
  setnames(graph, c("method", "value", "n"), c("method", "value", "n"))
  graph[, method := unname(labels[method])]
  graph[, percent := 100 * n / sum(n), by = method]
  graph[, method := factor(method, levels = rev(display_order))]
  fwrite(graph, file.path(output_dir, paste0(name, "_halfeye_data.csv")))
  variable <- graph[, .(k = uniqueN(value)), by = method][k > 1, method]
  # A degenerate distribution has no density; retain its zero-width boxplot.
  p <- ggplot2::ggplot(graph, ggplot2::aes(x = method, y = value, fill = method, weight = percent)) +
    ggdist::stat_halfeye(data = graph[method %in% variable],
      adjust = if (name %in% c("new_sim2","new_sim2_2")) 5 else .5,
      justification = -.2, .width = 0, point_colour = NA, show.legend = FALSE) +
    ggplot2::geom_boxplot(width = .12, outlier.color = NA, alpha = .5, show.legend = FALSE) +
    ggplot2::scale_x_discrete(limits=rev(display_order)) +
    ggplot2::scale_fill_manual(values=setNames(scales::hue_pal()(5),sort(unname(labels))),drop=FALSE) +
    ggplot2::labs(x = "Method", y = "Raw Bias") +
    # Zoom after density estimation: do not discard tails outside the display window.
    ggplot2::coord_flip(ylim = c(-.5, .5)) +
    jtools::theme_apa()
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".pdf")), p,
    device = grDevices::pdf, width = 6, height = 6, useDingbats = FALSE)
  if (requireNamespace("ragg", quietly = TRUE))
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), p,
      device = ragg::agg_png, width = 6, height = 6, dpi = 300)
  if (save_tex) {
    tikzDevice::tikz(file.path(output_dir, paste0(name, ".tex")),
                     width = 5.5, height = 5, standAlone = TRUE)
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  invisible(p)
}



if (DRAW_PLOTS) for (sn in unique(frequency$scenario)) {
  graph1 <- frequency[scenario==sn,.(n=sum(n)),by=.(method=sub("^bias_","",estimator),value)]
  name <- switch(sn,general="new_sim1",common="new_sim2",ignorability="new_sim3")
  plot_halfeye(graph1,name,OUTPUT_DIR)
}

if (DRAW_PLOTS) {
  labels <- c(bias_naive="Naive",bias_ancova="ANCOVA",bias_change="Change (F)",bias_change_2="Change (A)",bias_change_1="Change (D)")
  g <- copy(curves)
  g[,method:=factor(unname(labels[estimator]),levels=unname(labels))]
  setorder(g,method,diff)
  endpoints <- g[, .SD[.N], by=method]
  setorder(endpoints,mean_abs_bias)
  endpoints[,label_y:=mean_abs_bias]
  gap <- max(base::diff(range(g$mean_abs_bias))*.07,.001)
  if(nrow(endpoints)>1) for(k in 2:nrow(endpoints)) endpoints$label_y[k] <- max(endpoints$label_y[k],endpoints$label_y[k-1]+gap)
  endpoints[,label_x:=diff+max(g$diff)*.06]
  p <- ggplot2::ggplot(g,ggplot2::aes(diff,mean_abs_bias,linetype=method))+
    ggplot2::geom_line(linewidth=.85)+
    ggplot2::geom_segment(data=endpoints,ggplot2::aes(x=diff,xend=label_x,y=mean_abs_bias,yend=label_y),inherit.aes=FALSE,linewidth=.3)+
    ggplot2::geom_text(data=endpoints,ggplot2::aes(x=label_x,y=label_y,label=method),inherit.aes=FALSE,hjust=0,size=3.3,fontface="italic")+
    ggplot2::scale_linetype_manual(values=c("solid","dashed","dotted","dotdash","longdash"))+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.02,.35)))+
    ggplot2::labs(x="Violation of Ignorability Assumption",y="Mean Absolute Asymptotic Bias")+
    jtools::theme_apa()+ggplot2::theme(legend.position="none")
  ggplot2::ggsave(file.path(OUTPUT_DIR,"sim_1_ig_2.pdf"),p,device=grDevices::pdf,width=6,height=6,useDingbats=FALSE)
  if(requireNamespace("ragg",quietly=TRUE)) ggplot2::ggsave(file.path(OUTPUT_DIR,"sim_1_ig_2.png"),p,device=ragg::agg_png,width=6,height=6,dpi=300)
}

writeLines(c(paste("Grid:",paste(ng,collapse=", ")),paste("Inner dimensions:",INNER_DIM),
  "Complete enumeration of free coefficients; common-trend derived coefficients may be off-grid.",
  capture.output(sessionInfo())),file.path(OUTPUT_DIR,"session_info.txt"))
print(summary_table)
message("Saved to: ",normalizePath(OUTPUT_DIR))
