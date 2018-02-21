# Graphing vertical migration rates of genus' from `Feeley et al. 2011` (DOI: 10.1111/j.1365-2699.2010.02444.x)

# Packages
library(ggplot2)
library(dplyr)

# Read in Migration rate (m yr^-1)
Migration_rate <- c(5.4,2.4,5.7,8.1,7.7,-1.9,-0.2,5.4,-0.1,3.5,24.9,9.4,-0.3,
                    0.2,-0.5,1.7,5.5,0.8,-1.3,-6.2,-1.0,0.2,0.8,-1.3,0.0,3.7,
                    1.7,3.4,-5.3,-1.1,-1.3,30.0,2.2,0.5,-0.2,0.6,0.3,-4.4)
## Check length
length(Migration_rate)

## Order data
Migration_rate_sorted <- sort(Migration_rate)

# Create Rank id vector
Rank <- seq(1,38, by = 1)

# Create data frame
Migration <- data.frame(Migration_rate_sorted, Rank)
head(Migration)

# Plot bar graph
ggplot(data = Migration, aes(x = Rank, y = Migration_rate_sorted)) + 
  geom_bar(stat = "identity", fill = "black", colour = "black", width = 0.75) +
  geom_hline(aes(yintercept=0), linetype = 5) +
  xlab("Genus") +
  ylab(expression(paste("Migration Rate (m yr"^-1,")"))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  annotate("text", x = 9, y = -4, label = " = 14", size = 15) +
  annotate("text", x = 30, y = 20, label = " = 22", size = 15) + 
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25))
  

