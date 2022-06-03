#' Calculates isolate group membership
#'
#' This function calculates the euclidean distance from ancestry matrices (Q). It assigns isolates to groups using a pairwise distance of 0 (+/- 1 SD). These membership groups are then plotted (if plot = TRUE) as a network and output as a membership table at each specified K value.
#' @param K The K value (or range) used to calculate ancestry matrix
#' @param plot TRUE/FALSE: Plots relationships between isolates visually
#' @return membership.txt file containing membership assignments at each K value
#' @return Kplot.tiff TIFF file containing network plot at specified K value
#' @export
#' @examples
#' calc_dist(K = 4)
#' calc_dist(K = 2:20, plot = FALSE)

calc_dist <- function(K, plot = TRUE){

minK = K[1]

#test minimum K greater than 1
if (minK < 2) {
  stop("'K' can't be less than 2")
}

#reads in data from file list
making_tbls <- function(file_list) {
  my_data <- list()
  for (i in seq_along(file_list)) {
    my_data[[i]] <- read.table(file = paste("q_files/", file_list[i], sep = ""), header = F)
  }
   return(my_data)
}


labels <- read.table("label.txt") # reads label data from VCF file

# wrangles data into correct format
calc_euclid <- function(data_list) {
  input <- cbind(labels[1], data_list)
  euclid <- suppressWarnings(dist(input, method = "euclidean"))
  euclid <- as.matrix(euclid, labels = T)
  colnames(euclid) = rownames(euclid) = input$V1
  euclid_long <- reshape2::melt(euclid)
  return(euclid_long)
}

#create empty lists
file.list = list()
euclid = list()
merged_data = list()
merged_data_stats = list()
k_membership = list()

# calculating euclidean distances, plotting networks and outputting membership assignment
for (k in K) {
  print("*************************************")
  p = paste("Calculating Euclidean distance for K =", k)
  print(p)
  print("*************************************")

  # make file lists
  file.list[[k]] <- list.files(path = "q_files", pattern = paste0(k,".Q"), full.names = FALSE)
  # create list of data
  data <- lapply(file.list, making_tbls)
  #calculate euclidean distances from data
  euclid[[k]] <- lapply(data[[k]], FUN = calc_euclid)
  #Merge all dataframes of within list into a single dataframe
  merged_data[[k]] = suppressWarnings(Reduce(function(...) merge(..., all=T, by=c(1,2)), euclid[[k]]))
  colnames(merged_data[[k]]) = c("Var1", "Var2", 3:length(merged_data[[k]]))
  #creating new dataframe containing statistics
  merged_data_stats[[k]] <- merged_data[[k]][1:2]
  #calculate the mean
  merged_data_stats[[k]]$mean = rowSums(merged_data[[k]][3:length(merged_data[[k]])]) / (length(merged_data[[k]])-2)
  #calculating the standard deviation
  merged_data_stats[[k]]$SD =  apply(merged_data[[k]][3:length(merged_data[[k]])], MARGIN = 1, FUN = sd)
  # empty unused K elements
  for (l in (minK-1):1) {
    merged_data_stats[[l]] = list()
  }

  #Plotting network analysis
  K_grouping <- subset(merged_data_stats[[k]], mean - 2*SD <= 0) #grouping based on 2 SD away from 0. Should be only 5% false negative, according to 68:95:99 rule.

  #Export image of network
  if (plot == TRUE) {
    nodelist <- unique.data.frame(as.data.frame(K_grouping[, 1]))
    nodelist$Size <- 2
    edgelist <- K_grouping[, c(1:3)]
    g <- igraph::graph_from_data_frame(edgelist, vertices = nodelist, directed = FALSE)
    #plot
    print("Exporting TIFF image")
    tiff(file=paste0("K", k,"plot.tiff"),width = 2000, height = 2000,res=300)
    plt <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link0(width=0.08,colour="black") +
      ggraph::geom_node_point(col="Black",size=1, position = 'identity') +
      ggraph::geom_node_text(aes(label=name), size = 3, color='black', repel = TRUE) +
      ggraph::theme_graph()
    print(plt)
    dev.off()
  }

  #creating dataframe with all membership groups
  K_memb <- igraph::graph.data.frame(K_grouping)
  K_memb <- as.data.frame(igraph::clusters(K_memb)$membership)
  colnames(K_memb) <- paste0("Group K", k)

  # membership list
  k_membership[[k]] <- K_memb

  }
  final <- dplyr::bind_cols(k_membership)
  print("Exporting membership assignments")
  # export membership assignments
  write.table(final, file = paste0("membership.txt"), sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)

  ## creating recommendation for K value based on minimum standard deviation
  sd_list <- list()
  for (k in 2:20){
    sd_list[k] <- mean(merged_data_stats[[k]]$SD)
  }
  recommended_k <- which.min(unlist(sd_list))
  recommended_k_isolates <- final[recommended_k]
  write.table(recommended_k_isolates, file = paste0("membership_recommended.txt"), sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)
}


