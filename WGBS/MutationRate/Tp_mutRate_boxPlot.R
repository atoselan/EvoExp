# Make boxplot for Tp mutation rate

library(ggplot2)
library(RColorBrewer)


group <- c("T70_22C", "T70_22C", "T70_22C", "T70_22C", "T70_22C",
           "T270_22C", "T270_22C", "T270_22C", "T270_22C",
           "T32_9C", "T32_9C", "T32_9C", "T32_9C", "T32_9C",
           "T144_9C", "T144_9C", "T144_9C", "T144_9C", "T144_9C",
           "T250_9C", "T250_9C", "T250_9C", "T250_9C", "T250_9C",
           "T50_32C", "T50_32C", "T50_32C", "T50_32C", "T50_32C",
           "T210_32C", "T210_32C", "T210_32C", "T210_32C", "T210_32C",
           "T300_32C",
           "T450_32C", "T450_32C", "T450_32C")

sample <- c("T70_22C_TP1", "T70_22C_TP2", "T70_22C_TP3", "T70_22C_TP4", "T70_22C_TP5",
            "T270_22C_TP2", "T270_22C_TP3", "T270_22C_TP4", "T270_22C_TP5",
            "T32_9C_TP1", "T32_9C_TP2", "T32_9C_TP3", "T32_9C_TP4", "T32_9C_TP5",
            "T144_9C_TP1", "T144_9C_TP2", "T144_9C_TP3", "T144_9C_TP4", "T144_9C_TP5",
            "T250_9C_TP1", "T250_9C_TP2", "T250_9C_TP3", "T250_9C_TP4", "T250_9C_TP5",
            "T50_32C_TP1", "T50_32C_TP2", "T50_32C_TP3", "T50_32C_TP4", "T50_32C_TP5",
            "T210_32C_TP1", "T210_32C_TP2", "T210_32C_TP3", "T210_32C_TP4", "T210_32C_TP5",
            "T300_32C_TP1", 
            "T450_32C_TP2", "T450_32C_TP3", "T450_32C_TP5")

mutRate <- c(0, 1.12931E-09, 0, 0, 1.86005E-09,
             5.64657E-10, 5.45186E-10, 0, 0,
             1.37482E-09, 0, 0, 0, 5.8557E-10,
             8.54616E-10, 2.08032E-09, 8.1079E-10, 3.9526E-10, 3.85619E-10,
             1.50575E-09, 3.30915E-09, 2.51529E-09, 1.75671E-09, 6.87409E-10,
             0, 0, 7.9052E-10, 7.52876E-10, 2.87462E-09,
             9.8815E-10, 1.43731E-09, 4.18511E-09, 2.25863E-09, 0,
             2.63507E-09, 3.29383E-09, 8.38919E-09, 5.89015E-09)

# Combine to make data frame
df <- data.frame(group, sample, mutRate)


# Order groups for plotting
df$order <- ordered(df$group,
                    levels = c("T70_22C","T270_22C",
                               "T32_9C", "T144_9C", "T250_9C",
                               "T50_32C", "T210_32C", "T300_32C", "T450_32C"))


# Make colour palette
greenScale <- brewer.pal(9,"Greens")
blueScale <- brewer.pal(9,"Blues")
redScale <- brewer.pal(9,"Reds")

cols <- c(greenScale[4],
          greenScale[6],
          blueScale[2],
          blueScale[4],
          blueScale[6],
          redScale[2],
          redScale[4],
          redScale[6],
          redScale[8])


# Make boxplot
pdf("D://BCF_VAR_CALLS/noMisMatchFilteringAllSamples/MutationRate/Tp_mutRate.pdf")
  ggplot(df, aes(x=order, y=mutRate)) + 
    theme_bw() +
    geom_boxplot(outlier.shape=NA, notch=FALSE, fill=cols) +
    #geom_jitter(width=0.1, size=0.5, alpha=0.5) +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14, angle = 90),
          axis.text.y = element_text(size = 14),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    ylab("Mutation rate")
dev.off()