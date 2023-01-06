library(data.table)
library(MicrobiomeGS2)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(particles)

currency_mets <- c("cpd00001","cpd00067","cpd00002","cpd00003","cpd00004",
                   "cpd00009","cpd00008","cpd00006","cpd00005","cpd00011",
                   "cpd00010","cpd00013","cpd00012","cpd00023","cpd00013",
                   "cpd00007","cpd00014")

modfiles <- dir("data2/", full.names = TRUE, pattern = "\\.RDS$")

models <- lapply(modfiles, readRDS)

# rm dead ends
models <- lapply(models, function(x) {
  der <- deadEndMetabolites(x)
  x <- rmReact(x, react = der$der)
  return(x)
})



modj <- join_mult_models(models, merge.lb.method = "median")

# build network
nodes <- data.table(id = modj$modj@met_id,
                    name = modj$modj@met_name)
nodes[grepl("^M", id), comp := str_extract(id,"M[0-9]")]
nodes[is.na(comp), comp := "EX"]

edges <- list()

n <- dim(modj$modj@S)[2]

for(i in 1:n) {
  cat("\r",i,"/",n)
  rxn <- modj$modj@react_id[i]
  
  mts_LHS <- modj$modj@met_id[which(modj$modj@S[,i] < 0)]
  mts_RHS <- modj$modj@met_id[which(modj$modj@S[,i] > 0)]
  
  if(length(c(mts_LHS, mts_RHS)) > 1 & length(c(mts_LHS, mts_RHS)) < 12) {
    combis <- expand.grid(mts_LHS,mts_RHS)
    
    edges[[rxn]] <- data.table(from = combis$Var1,
                               to = combis$Var2,
                               rxn = rxn)
  }

  
}
edges <- rbindlist(edges)
edges <- edges[!grepl(paste(currency_mets, collapse = "|"), from)]
edges <- edges[!grepl(paste(currency_mets, collapse = "|"), to)]
edges[, compE := str_extract(rxn,"M[0-9]")]
edges[grepl("EX", rxn), compE := "EX"]

nodes <- nodes[id %in% c(edges$from, edges$to)]


#sort(table(gsub("M[0-9]_","",c(edges$from, edges$to))), decreasing = TRUE)

# makeing the network
nodes[, id := gsub("\\[|\\]","",id)]
edges[, to := gsub("\\[|\\]","",to)]
edges[, from := gsub("\\[|\\]","",from)]
netdata <- list(nodes = nodes,
                links = edges)
netdata <- as_tbl_graph(netdata, directed = FALSE, node_key = "id")


savings <- list()
k <- 1

save_status <- function(sim) {
  gr <- as_tbl_graph(sim)
  savings[[k]] <<- gr
  cat("\r",k)
  k <<- k + 1
}


simulation <- netdata |> 
  simulate() |> 
  wield(link_force) |> 
  wield(manybody_force) |> 
  wield(center_force) |> 
  wield(collision_force, radius = 0.5, n_iter = 2) |> 
  evolve(step = 255, on_generation = save_status)


# graph <- netdata |> 
#   simulate() |> 
#   wield(link_force) |> 
#   wield(manybody_force) |> 
#   wield(center_force) |> 
#   wield(collision_force, radius = 0.5, n_iter = 2) |> 
#   evolve() |> 
#   as_tbl_graph()

ncolors <- c(EX = "#FFFFFF",
             M1 = "#648FFF",
             M2 = "#785EF0",
             M3 = "#DC267F",
             M4 = "#FE6100",
             M5 = "#FFB000")


for(i in 1:length(savings)) {
  p <- ggraph(savings[[i]], 'nicely') + 
    geom_edge_link(aes(colour = compE), alpha = 0.2) + 
    geom_node_point(aes(fill = comp, color = comp), shape = 21, size = 1, alpha = 0.5,
                    stroke = 0.5) + 
    scale_fill_manual(values = ncolors) + 
    scale_color_manual(values = ncolors) +
    scale_edge_color_manual(values = ncolors) +
    #scale_edge_width('Value', range = c(0.5, 3)) + 
    coord_cartesian(ylim = c(-3000, 3192), xlim = c(-8000,3008), expand = FALSE) +
    theme_graph() +
    theme(legend.position = "none",
          aspect.ratio = 9/16,
          panel.background = element_rect(fill = "#0e0e0e", linewidth = 0),
          plot.background = element_rect(fill = "#0e0e0e", linewidth = 0),
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.spacing=grid::unit(c(0,0,0,0), "mm"))
  
  i_long <- as.character(i)
  i_long <- paste0(paste(rep("0", floor(log10(length(savings)))-floor(log10(i))),
                         collapse = ""),i_long)
  
  cat("\r",i_long)
  ggsave(paste0("frames/out_",i_long,".png"), plot = p, width =  16, height = 9,
         units = "in")
  if(i == length(savings)) {
    ggsave(paste0("frames/out_",paste0(rep("0",floor(log10(length(savings)))+1), collapse = ""),".png"), plot = p, width =  16, height = 9,
           units = "in")
  }
    
}


# ffmpeg -r 25 -f image2 -s 1920x1080 -i out_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ../test.mp4

p <- ggraph(savings[[255]], 'nicely') + 
  geom_edge_link(aes(colour = compE), alpha = 0.2) + 
  geom_node_point(aes(fill = comp, color = comp), shape = 21, size = 1, alpha = 0.5,
                  stroke = 0.5) + 
  scale_fill_manual(values = ncolors) + 
  scale_color_manual(values = ncolors) +
  scale_edge_color_manual('comp', values = ncolors) +
  #scale_edge_width('Value', range = c(0.5, 3)) + 
  coord_fixed(ylim = c(-3000, 3000), xlim = c(-8000,3000)) +
  theme_graph() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#0e0e0e"))
p


