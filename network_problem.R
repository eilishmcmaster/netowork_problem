library(sf)
library(sfnetworks)
library(TSP)
library(sf)
library(tidygraph)
library(igraph)
library(dplyr)
library(tibble)
library(ggplot2)
library(units)
library(ggmap)
library(openxlsx)
library(plotly)

#### import data ####
sp <- data.frame(X = 142.2004452133814, Y = -34.17478224432166, recordID='MILDURA')
starting_point <- st_as_sf(sp, coords = c("X", "Y"), crs = 4326) %>% st_transform(.,crs='EPSG:3577')

points <- read.xlsx('data/target_points.xlsx')
points_sf <- st_as_sf(points, coords = c("decimalLongitude", "decimalLatitude"), crs = 'EPSG:4326') %>% st_transform(crs = 'EPSG:3577')

points_sf <- bind_rows(points_sf, starting_point) 
#### site clusters ####
expand_m <- 100

# convex hulls for each cluster
hull_sf <- points_sf[!is.na(points_sf$cluster), ] %>%
  group_by(cluster) %>%
  slice(chull(st_coordinates(geometry))) %>%
  st_as_sf(coords = c("X", "Y"), crs = 'EPSG:3577') %>%  
  summarise(geometry = st_combine(geometry)) %>%  
  st_convex_hull() %>%                           
  filter(!st_geometry_type(geometry) %in% c("POINT")) %>%  
  st_buffer(dist = expand_m) %>%                 
  mutate(area = st_area(geometry))        

#### get roads ####
roads2 <- sf::read_sf('data/edited_cropped_roads_nonnull_valid.gpkg')
any(st_is_empty(roads2))
roads_linestring <- st_cast(st_cast(roads2, "MULTILINESTRING"),"LINESTRING")
roads_network <- as_sfnetwork(roads_linestring,directed = FALSE) %>% st_transform(.,crs='EPSG:3577')

#### find the closest point in each cluster (site) to the roads ####
distances <- st_distance(points_sf, roads_network)
min_distances <- apply(distances, 1, min)
min_dist <- data.frame(recordID = points_sf$recordID, min_distance = min_distances,
                       cluster = points_sf$cluster)

min_dist2 <- min_dist %>%
  filter(!is.na(cluster)) %>% 
  group_by(cluster) %>%
  filter(min_distance == min(min_distance, na.rm = TRUE)) %>%  # Find minimum
  ungroup() %>%
  bind_rows(min_dist %>% filter(is.na(cluster))) 

closest_points <- points_sf[points_sf$recordID %in% min_dist2$recordID,]

#### blend nodes into road network ####
#put the target sites on the road as nodes
blended <- st_network_blend(x=roads_network, closest_points, tolerance = 3000) # change tolerance to flter points far from roads
blended_nodes <- st_as_sf(activate(blended,'nodes'))
target_nodes <- blended_nodes[!is.na(blended_nodes$recordID),]
targeted_records <- points_sf[(points_sf$recordID %in% target_nodes$recordID),]
blended_edges <- st_as_sf(activate(blended,'edges'))

# Define the weight system for each road class
road_weights <- c("Principal Road" = 1, "Secondary Road" = 1.25, "Minor Road" = 1.5, "Track" = 2, "Track2" = 2)
blended_edges$weight <- road_weights[blended_edges$CLASS]*st_length(blended_edges)
blended <- activate(blended, 'edges') %>%  mutate(weight = blended_edges$weight)


####### HERE IS THE PROBLEM, TWO POINTS ARE UNREACHABLE ##########
# calculate the distance matrix using the updated weights
distance_matrix <- st_network_cost(blended, 
                                   from = target_nodes, to = target_nodes,
                                   # weights = blended$weight,
                                   weights=NULL,
                                   Inf_as_NaN = TRUE, direction = 'all')

distance_matrix2 <- units::drop_units(distance_matrix%>% as.matrix())

# workaround to drop the unreachable points from the TSP
index_reachable <- which(colSums(distance_matrix2, na.rm=TRUE)!=0)
distance_matrix3<- distance_matrix2[index_reachable,index_reachable]
tsp_problem <- TSP(distance_matrix3)

unreachable_records <- target_nodes[which(colSums(distance_matrix2, na.rm=TRUE)==0),]

# solve the TSP
tour <- solve_TSP(tsp_problem, method = "cheapest_insertion", weight=NULL, 
                  start = which(target_nodes$recordID=='MILDURA'))
optimal_path <- labels(tour)
ordered_points <- target_nodes[optimal_path, ]

# loop through consecutive points and get the route between them
routes <- list()
for(i in 1:(nrow(ordered_points) - 1)) {
  route <- st_network_paths(
    blended,
    weights = blended$weight,
    from = ordered_points[i,],
    to = ordered_points[i + 1,]
  )
  routes[[i]] <- route$edge_paths[[1]]
}

# extract coordinates and convert to data frames
target_nodes_coords <- st_coordinates(target_nodes) # Get coordinates for target_nodes
targeted_records_coords <- st_coordinates(targeted_records) # Get coordinates for targeted_records

# Convert to data frames and add recordID and coordinates
target_nodes_df <- as.data.frame(target_nodes) %>%
  mutate(X_start = target_nodes_coords[, "X"],
         Y_start = target_nodes_coords[, "Y"])

targeted_records_df <- as.data.frame(targeted_records) %>%
  mutate(X_end = targeted_records_coords[, "X"],
         Y_end = targeted_records_coords[, "Y"])

# merge the data frames by recordID
lines_df <- merge(target_nodes_df, targeted_records_df, by = "recordID")


# plot shows the route calculated. black points are target points for routing. the two red points are the unreachable points. 
x <- ggplot() +
  geom_sf(data = hull_sf, fill='tan2',color='tan2')+
  geom_sf(data = blended_edges, color = 'grey90') +
  geom_sf(data = blended_edges[unlist(routes), ], mapping = aes(color = CLASS)) + 
  geom_sf(data = target_nodes) +
  geom_sf(data = targeted_records, color = 'steelblue') +
  geom_sf(data = unreachable_records, color = 'red') +
  geom_segment(data = lines_df,
               aes(x = X_start, y = Y_start, xend = X_end, yend = Y_end),
               color = "black") +
  geom_sf(data = ordered_points[c(1, nrow(ordered_points)), ], color = "green", size = 2) +
  coord_sf(crs=3577, datum=3577)+
  theme_void()

plotly::ggplotly(x)
