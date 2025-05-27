# R coding basic skills -----------------------------------------------------
# Written for the Cal-Neva AFS workshop 2025-05-30 in Lodi
# Author: Eric Holmes, contact: eric.holmes@water.ca.gov

# Skills for new coders: 
#         1) learn how to trouble shoot your code: 
#             a. reousrces: help documents, forums, AI chatbots
#             b. ask the right questions: be explicit
#             c. know how to make a minimal reproducible example
#         2) package knowledge
#         3) Break a large project into components and work step by step:
#             a. create an outline to write code like you would write an essay

## Load Libraries -----------------------------------------------------------

library(tidyverse)
library(sf)
library(leaflet)

##  Indexing ----------------------------------------------------------------

df <- data.frame("col1" = letters, "col2" = 1:26)

#First row and first column
df[1,1]
#First row and all columns
df[1,]
#All rows and second column
df[,2]

## Operators ----------------------------------------------------------------
# == != > < >= | & %in% <- + - ^

##Equal to
df[df$col1 == "a",]
##Not equal to
df[df$col1 != "a",]
##Less than or equal to
df[df$col2 <= 5,]
##Greater than or equal to
df[df$col2 >= 5,]
## | Element-wise OR
df[df$col1 == "a" | df$col2 == 2,]
## & Element-wise AND
df[df$col1 == "a" & df$col2 == 1,]
##Membership
df[df$col1 %in% c(letters[1:5]),]

## Loops --------------------------------------------------------------------

##For loop
num_vector <- df$col2
for(i in num_vector){
  print(i)
}

##While loop
count = 1
while(count <= 26){
  print(count)
  count <- count + 1
}

##Using a loop, try to build a dataframe identical to our original df
df2 <- data.frame()
for(j in 1:26){
  print(j)
  df2 <- rbind(df2, data.frame(col1 = letters[j], col2 = j))
}

##Check if it worked
df == df2

## Functions ----------------------------------------------------------------

## Make a function which creates a dataframe identical to our original df
create_df <- function(rownums = 1){
  return(data.frame(col1 = letters[rownums],
                    col2 = rownums))
}

##run the function
df3 = create_df(1:26)

##Check if it worked
df == df3

## Conditionals -------------------------------------------------------------

df$col3 <- ifelse(df$col1 %in% letters[1:13], "strawberry", "banana")

## nested 
df$col4 <- ifelse(df$col1 %in% letters[1:8], "strawberry", 
                  ifelse(df$col1 %in% letters[9:16], "banana", 
                         ifelse(df$col1 %in% letters[17:24], "smoothie", "RACCOONS")))

##Using a conditional to see if df4 exists. If so, delete it; if not, create it
if(exists("df4")){print("df4 already exists, removing df4")
  rm(df4)} else{
    print("df4 cannot be found, creating df4")
    df4 = create_df(1:26)}

## dplyr: Split-apply-combine -----------------------------------------------

df5 <- data.frame(col1 = rep(letters[1:5], each = 100),
                  col2 = rep(1:5, each = 100) + rnorm(n = 5*100, mean = 1, sd = 1))

ggplot(df5, aes(x = col1, y = col2)) +  
  geom_boxplot(fill = "cyan3") + geom_jitter(height = 0, shape = 1) + theme_bw()

df5ply <- df5 %>% group_by(col1) %>% 
  summarize(avg2 = mean(col2), stdev2 = sd(col2))

ggplot(df5, aes(x = col1)) +  
  geom_jitter(aes(y = col2), height = 0, shape = 1, alpha = .2) + theme_bw() +
  geom_point(data = df5ply, aes(y = avg2)) + 
  geom_errorbar(data = df5ply, aes(ymin = avg2 - stdev2, ymax = avg2 + stdev2), width = .1)

## tidyr: reshape data ------------------------------------------------------

df5$ID <- rep(paste0("xs", 1:100), 5)

df5wide <- df5 %>% pivot_wider(names_from = "col1", values_from = "col2")

df5long <- df5wide %>% pivot_longer(cols = c(letters[1:5]),
                                    names_to = "col1", values_to = "col2")

## Joining/merging data -----------------------------------------------------

key <- data.frame(col1 = letters[1:5], fruit = c("apple", "banana", "cantelope", "date", "elderberry"))

#left merge
df5merge <- merge(df5, key, by = "col1", all.x = T)

## Basic mapping tutorial with simple features (sf) package -----------------

# Define domain for coordinates
xmin <- -122
xmax <- -121
ymin <- 38
ymax <- 39

#Create 100 random points within domain
rand_pts <- data.frame(id = 1:100,
                       lon = runif(100, xmin, xmax),
                       lat = runif(100, ymin, ymax))

# Create an sf point object
rand_pts_sf <- st_as_sf(rand_pts, coords = c("lon", "lat"), crs = 4326)

# Create a matrix of coordinates for a box polygon
boxdf <- data.frame(lon = c(xmin, xmax, xmax, xmin, xmin),
                    lat = c(ymin, ymin, ymax, ymax, ymin))

# Create an sf polygon object
box_polygon <- st_polygon(list(as.matrix(boxdf)))

# Convert polygon to an sf object
box_sf <- st_sf(geometry = st_sfc(box_polygon), crs = 4326)

## Simple plot
plot(box_sf$geometry)
plot(rand_pts_sf$geometry, add = T)

## simple ggplot
ggplot() + 
  geom_polygon(data = boxdf, aes(x = lon, y = lat), fill = NA, color = "red") + 
  geom_point(data = rand_pts, aes(x = lon, y = lat)) + 
  theme_minimal()

## simple ggplot with projected coordinate space
ggplot() + geom_sf(data = box_sf, fill = NA) + 
  geom_sf(data = rand_pts_sf) + theme_minimal()

## Plot sampling sites in an interactive map
leaflet() %>% addPolygons(data = box_sf, fill = NA) %>% 
  addCircles(data = rand_pts_sf) %>% addTiles()
