#### Colorful R Plots with Wes Anderson Palettes & ggplot2 - Pirate Ship Features
# data.world/funsizemaddy/pirate2/workspace/file?filename=Pirates+%281%29.xlsx

# set the working directory
#setwd("/Users/bamflappy/Documents/rPlayground")

# import the pirate ship data from the PiratesShip.csv file
ships <- read.csv("PiratesShip.csv")

# view the column names of the ship data
colnames(ships)

# view the first few lines of ship data
head(ships)

# explore the relationship between Sails and CrewCapacity
# first, create a scatter plot
ggplot(data = ships, aes(x = Sails, y = CrewCapacity)) +
  geom_point()

ggsave("ship_plot_crew_sails.png", plot = last_plot())

# color the scatter plot by MaidenYear
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point()

ggsave("ship_plot_crew_sails_year.png", plot = last_plot())

# create a colorful plot by searching the internet for "ggplot wes anderson"
# https://github.com/karthik/wesanderson
# https://rforpoliticalscience.com/2020/07/26/make-wes-anderson-themed-graphs-with-wesanderson-package-in-r/

# if not already, install the wesanderson package
#install.packages("wesanderson")

# load the wesanderson package
library(wesanderson)

# checkout the color palette options
names(wes_palettes)

# look at the info for the wes_palette function
?wes_palette

# look at the info for the scale_color_gradientn function
?scale_color_gradientn

# color the scatter plot using wes_palettes and scale_color_gradientn
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))

ggsave("ship_plot_crew_sails_year_color.png", plot = last_plot())

# add axis and plot titles
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title="Pirate Ship Crew Capacity by Sails", 
       x ="Number of Sails", 
       y = "Crew Capacity")

ggsave("ship_plot_crew_sails_year_title.png", plot = last_plot())

# adjust the axis and plot title colors
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title="Pirate Ship Crew Capacity by Sails", 
       x ="Number of Sails", 
       y = "Crew Capacity") +
  theme(
    plot.title = element_text(color="darkgreen", size=14, face="bold.italic"),
    axis.title.x = element_text(color="darkblue", size=14, face="bold"),
    axis.title.y = element_text(color="darkred", size=14, face="bold")
  )

ggsave("ship_plot_crew_sails_year_colorTitle.png", plot = last_plot())

# center the plot title
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title="Pirate Ship Crew Capacity by Sails", 
       x ="Number of Sails", 
       y = "Crew Capacity") +
  theme(
    plot.title = element_text(color="darkgreen", size=14, face="bold.italic"),
    axis.title.x = element_text(color="darkblue", size=14, face="bold"),
    axis.title.y = element_text(color="darkred", size=14, face="bold")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# save the last plot using the ggsave function
ggsave("ship_plot_crew_sails_year_centerTitle.png", plot = last_plot())
