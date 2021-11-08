# Script to generate individual figures for all scenarios
# 
# 
# David J Price
# 22 October, 2021



outer.alpha <- 0.15
inner.alpha <- 0.4
line.alpha <- 0.8

total.inf.upper.y.nsw <- 250
arrivals.inf.upper.y.nsw <- 250

total.inf.upper.y.nsw.baseline <- 25000
arrivals.inf.upper.y.nsw.baseline <- 25000

total.inf.upper.y.wa <- 5000
arrivals.inf.upper.y.wa <- 5000


abm.folder <- "~/Dropbox/wp3outputs/"
all.files <- list.files(abm.folder)

all.files <- all.files[grepl(pattern = ".csv", x = all.files)]

for (the.file in all.files){
  
  dat <- data.table::fread(paste0(abm.folder, the.file))
  

  dat <- dat %>% mutate(Q = case_when(
    Q == 0 ~ 0,
    Q > 0 ~ 1)
  ) %>% complete(Date = 1:500, Q = 0:1, Sim = 1:200, fill = list(n = 0))
  
  plot.dat.total <- dat %>% 
    group_by(Date, Sim) %>% 
    mutate(n = sum(n)) %>% 
    group_by(Date) %>% 
    summarise(x=list(enframe(quantile(n, probs=c(0.05, 0.25, 0.75, 0.95)), "quantiles", "value"))) %>%
    unnest(x) %>% mutate(quantiles = as.numeric(gsub(pattern = "%",replacement = "",x = quantiles))) %>%
    mutate(quantiles = paste0("q",quantiles)) %>% pivot_wider(names_from = quantiles, values_from=value)
  
  
  dat %>% ggplot() + aes(x = Date, y = n, group = Sim, colour = Q) + geom_line(alpha = 0.2)

  
  total.inf.plot <- plot.dat.total %>% 
    ggplot() + 
    aes(x = Date) + 
    geom_ribbon(aes(ymin = q5, ymax = q95), fill = "steelblue3", alpha = outer.alpha) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue3",alpha = inner.alpha) +
    scale_y_continuous("Daily New Cases") +
    scale_x_continuous("Day") +
    cowplot::theme_cowplot() +
    theme(legend.position = "top",
          text = element_text(size = 16))
  
  if(grepl("NSW", x = the.file)){
    if(grepl(pattern = "_baseline", x = the.file)){
      total.inf.plot <- total.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                         ylim = c(0, total.inf.upper.y.nsw.baseline))  
    } else{
      total.inf.plot <- total.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                         ylim = c(0, total.inf.upper.y.nsw))
    }
  } else{
    total.inf.plot <- total.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                       ylim = c(0, total.inf.upper.y.wa))
  }
  
  
  plot.dat.split <- dat %>% 
    mutate(`Cluster Index` = case_when(
      Q == 0 ~ "Local",
      Q > 0 ~ "Arrival"
    )) %>% 
    group_by(`Cluster Index`, Date, Sim) %>% 
    mutate(n = sum(n)) %>% 
    group_by(`Cluster Index`, Date) %>% 
    summarise(x=list(enframe(quantile(n, probs=c(0.05, 0.25, 0.75, 0.95)), "quantiles", "value"))) %>%
    unnest(x) %>% mutate(quantiles = as.numeric(gsub(pattern = "%",replacement = "",x = quantiles))) %>%
    mutate(quantiles = paste0("q",quantiles)) %>% pivot_wider(names_from = quantiles, values_from=value)
  
  arr.inf.plot <- plot.dat.split %>% 
    ggplot() + 
    aes(x = Date, fill = `Cluster Index`) + 
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha) +
    scale_y_continuous("Daily New Cases") +
    scale_x_continuous("Day") +
    # coord_cartesian(xlim = c(0, 365), ylim = c(0, upper.y*1.05)) +
    cowplot::theme_cowplot() +
    theme(legend.position = "top",
          text = element_text(size = 16))
  
  if(grepl("NSW", x = the.file)){
    if(grepl(pattern = "_baseline", x = the.file)){
      arr.inf.plot <- arr.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                     ylim = c(0, arrivals.inf.upper.y.nsw.baseline))
    } else{
      arr.inf.plot <- arr.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                     ylim = c(0, arrivals.inf.upper.y.nsw))
    }
  } else{
    arr.inf.plot <- arr.inf.plot + coord_cartesian(xlim = c(0, 365), 
                                                   ylim = c(0, arrivals.inf.upper.y.wa))
  }
  
  output.folder <- paste0("./output_25Oct2021/", gsub(pattern = ".csv", replacement = "", the.file), "/")
  dir.create(output.folder, showWarnings = FALSE)
  save(total.inf.plot, arr.inf.plot, file = paste0(output.folder,gsub(pattern = ".csv", replacement = "", the.file),"_plots.Rdata"))
  ggsave(paste0(output.folder, "inf_by_source.png"), plot = arr.inf.plot, height = 6, width = 8, units = "in", bg = "white")
  ggsave(paste0(output.folder, "inf_total.png"), plot = total.inf.plot, height = 6, width = 8, units = "in", bg = "white")
  
}


