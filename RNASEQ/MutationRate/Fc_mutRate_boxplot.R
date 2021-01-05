# Make boxplot for Fc mutation rate

library(ggplot2)


group <- c("T66_3.5C", "T66_3.5C", "T66_3.5C",
           "T47_8C", "T47_8C", "T47_8C")

sample <- c("T66_3pt5C_FC1", "T66_3pt5C_FC2", "T66_3pt5C_FC3",
            "T47_8C_FC1", "T47_8C_FC2", "T47_8C_FC1")

mutRate <- c(1.59905E-09, 2.25748E-09, 3.10403E-09,
             2.9059E-09, 3.03799E-09, 5.01929E-09)

# Combine to make data frame
df <- data.frame(group, sample, mutRate)

# Make boxplot
pdf("D://BCF_VAR_CALLS/noMisMatchFilteringAllSamples/MutationRate/FC_mutRate.pdf")
  ggplot(df, aes(x=reorder(group,mutRate), y=mutRate, fill=group)) + 
    theme_bw() +
    geom_boxplot(outlier.shape=NA, notch=FALSE) +
    #geom_jitter(width=0.1, size=0.5, alpha=0.5) +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          #legend.title = element_text(size = 16),
          #legend.text = element_text(size = 14),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    ylab("Mutation rate")
dev.off()
