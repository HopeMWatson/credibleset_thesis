#Figures and Plots for Thesis - Overcoverage

###############################################################################################################################
##General_ordering_barplot_size
General_ordering_barplot_size <- ggplot()

##General_ordering_barplot_coverage


###############################################################################################################################
#OR

##OR Size 
OR_ordering_barplot_size <- ggplot(df, aes(x=OR, y=size, fill=Ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=size-lci, ymax=size+uci), width=0.2, 
                position=position_dodge(0.9))

OR_ordering_barplot_size+labs(title "Size in Varying Odds Ratios", x="Odds Ratio", y="Size")+
  theme_minimal()
theme(legend.title=element_blank())

print(OR_ordering_barplot_size)


##OR Coverage 
OR_ordering_barplot_coverage <- ggplot(df, aes(x=OR, y=coverage, fill=ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=coverage-lci, ymax=coverage+uci), width=0.2, 
                position=position_dodge(0.9))

OR_ordering_barplot_coverage+labs(title "Coverage in Varying Odds Ratios", x="Odds Ratio", y="Coverage")+
  theme_minimal()
  theme(legend.title=element_blank())
  
print(OR_ordering_barplot_coverage)


###############################################################################################################################
#N

##N Size
N_ordering_barplot_size <- ggplot(df, aes(x=n, y=size, fill=Ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=size-lci, ymax=size+uci), width=0.2, 
                position=position_dodge(0.9))

N_ordering_barplot_size+labs(title "Size in Varying Sample Sizes", x="Sample Size (n)", y="Size")+
  theme_minimal()
  theme(legend.title=element_blank())

print(N_ordering_barplot_size)


##N Coverage
N_ordering_barplot_coverage <- ggplot(df, aes(x=n, y=coverage, fill=Ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=coverage-lci, ymax=coverage+uci), width=0.2, 
                position=position_dodge(0.9))

N_ordering_barplot_coverage+labs(title "Coverage in Varying Sample Sizes", x="Sample coverage (n)", y="coverage")+
  theme_minimal()
  theme(legend.title=element_blank())

print(N_ordering_barplot_coverage)




##############################################################################################################################
#Thr

##Thr Size

Thr_ordering_barplot_size <- ggplot(df, aes(x=thr, y=size, fill=Ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=size-lci, ymax=size+uci), width=0.2, 
                position=position_dodge(0.9))

Thr_ordering_barplot_size+labs(title "Size in Varying Thresholds", x="Threshold", y="Size")+
  theme_minimal()
  theme(legend.title=element_blank())

print(Thr_ordering_barplot_size)

##Thr Coverage 
Thr_ordering_barplot_coverage <- ggplot(df, aes(x=thr, y=coverage, fill=Ordered)) + 
  geom_bar (stat="identity", color="black", Position(position_dodge())) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2, 
                position=position_dodge(0.9))

Thr_ordering_barplot_coverage+labs(title "Coverage in Varying Thresholds", x="Threshold", y="coverage")+
  theme_minimal()
  theme(legend.title=element_blank())

print(Thr_ordering_barplot_coverage)