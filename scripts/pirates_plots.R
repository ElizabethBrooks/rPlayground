#### Fun R Plots with ggplot2 - Pirate Features
# data.world/funsizemaddy/pirate2/workspace/file?filename=Pirates+%281%29.xlsx

# set the working directory
#setwd("/Users/bamflappy/Documents/rPlayground")

# if not already, install the ggplot2 package
#install.packages("ggplot2")

# load the ggplot2 library
library(ggplot2)

# import the pirate data from the PiratesPirate.csv file
pirates <- read.csv("data/PiratesPirate.csv")

# view the column names of the pirates data
colnames(pirates)

# view the first few lines of pirates data
head(pirates)

# check out the basic ggplot function and geoms
# https://ggplot2.tidyverse.org/reference/
# https://datacarpentry.org/r-socialsci/04-ggplot2/index.html

### One Dimension - Teeth or Limbs

# start with plotting only the Limbs dimension of the pirates data
ggplot(data = pirates, aes(x = Limbs)) +
  geom_bar()

# next, plot only the Teeth dimension of the pirates data
ggplot(data = pirates, aes(x = Teeth)) +
  geom_bar()

### Two Dimensions - Teeth & Limbs

# create a scatter plot comparing the two dimensions of Limbs vs Teeth
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_point()

# modify the scatter plot to include the count of pirates at each point
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_count()

# check the data for the outlier pirate with 2.5 Limbs
pirates[pirates$Limbs == "2.5", ]

# change to a boxplot to better describe the data at each point
ggplot(data = pirates, aes(x = Limbs, y = Teeth)) +
  geom_boxplot()

# check the info for the factor function
?factor

# use the factor function to make the different numbers of Limbs into categories
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot()

# create a more detailed plot using colors to add another dimension
# look up color options on the internet by searching "ggplot boxplot color"
# example 1 from r-graph-gallery.com/264-control-ggplot2-boxplot-colors.html
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="red", fill="orange", alpha=0.2)

# example 2 from r-graph-gallery.com/264-control-ggplot2-boxplot-colors.html
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth, fill = Origin)) +
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none")

### Three Dimensions - Teeth & Limbs & Origin

# create a set of boxplots with one for each Origin of pirates
# look up facet options on the internet by searching "ggplot boxplot facet"
# www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin)

# change the appearance of the individual plot titles
# www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))

# change the axis titles and add a plot title
# look up how to add titles on the internet by searching "ggplot title"
# top result www.sthda.com/english/wiki/ggplot2-title-main-axis-and-legend-titles
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) +
  labs(title="Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth")

# adjust the colors of the axis and plot titles
# www.sthda.com/english/wiki/ggplot2-title-main-axis-and-legend-titles
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) +
  labs(title="Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth") +
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="purple", size=14, face="bold")
  )

# center the plot title
# look up how to center the title on the internet by searching "ggplot center title"
# https://stackoverflow.com/questions/40675778/center-plot-title-in-ggplot2
ggplot(data = pirates, aes(x = factor(Limbs), y = Teeth)) +
  geom_boxplot(color="darkgreen", fill="orange", alpha=0.2) +
  facet_wrap(~ Origin) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) +
  labs(title="Comparison of Pirate Teeth and Limb Numbers by Origin", 
       x ="Number of Limbs", 
       y = "Number of Teeth") +
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="purple", size=14, face="bold")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# check out the info for the ggsave function
?ggsave

# check out the info for the last_plot function
?last_plot

# save the last plot using the ggsave function
ggsave("plots/pirates_plot_teeth_limbs_origin.png", plot = last_plot())
