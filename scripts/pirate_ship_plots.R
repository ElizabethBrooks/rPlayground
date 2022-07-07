#### Colorful R Plots with Wes Anderson Palettes & ggplot2 - Pirate Ship Features

# if not already, install the ggplot2 package
#install.packages("ggplot2")

# load the ggplot2 library
library(ggplot2)

# import the pirate ship data from the PiratesShip.csv file
ships <- read.csv("data/PiratesShip.csv")

# view the column names of the ship data
colnames(ships)

# view the first few lines of ship data
head(ships)

### One Dimension of Data with Vectors - CrewCapacity or Sails
# create a plot for each dimension of CrewCapacity or Sails

# start with plotting only the CrewCapacity dimension of the ship data
ggplot(data = pirates, aes(x = CrewCapacity)) +
  geom_bar()

# next, plot only the Sails dimension of the ship data
ggplot(data = pirates, aes(x = Sails)) +
  geom_bar()

### Two Dimensions of Data with Dataframes - CrewCapacity & Sails
# combine the two dimensions into one plot with geom_point

# first, create a scatter plot
ggplot(data = ships, aes(x = Sails, y = CrewCapacity)) +
  geom_point()

### Three Dimensions of Data with Colors - CrewCapacity & Sails & MaidenYear
# add a third dimension to your plot with color

# color the scatter plot by MaidenYear
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point()

### Fun Colors - Wes Anderson Palette
# create a fun colorful plot using the wesanderson color palette

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

### Saving Plots - ggsave

# check out the info for the ggsave function
?ggsave

# check out the info for the last_plot function
?last_plot

# save the last plot using the ggsave function
ggsave("plots/ship_plot_crew_sails_year_color.png", plot = last_plot())

### Bonus Exercises - Adjusting Plot Appearance

# add axis and plot titles
ggplot(data = ships, aes(x = Sails, y = CrewCapacity, color = MaidenYear)) +
  geom_point() +
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  labs(title="Pirate Ship Crew Capacity by Sails", 
       x ="Number of Sails", 
       y = "Crew Capacity")

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
ggsave("plots/ship_plot_crew_sails_year_centerTitle.png", plot = last_plot())
