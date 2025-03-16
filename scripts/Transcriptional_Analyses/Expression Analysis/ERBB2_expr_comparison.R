# Compare ERBB2 expression between NaiveLike and TreatedLike?
setwd("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis")

load("Processed_ExprMatrix.Rdata")

# Extract TYMS expression
samples <- expr.matrix[expr.matrix$status == "NaiveLike.Sig17+" | expr.matrix$status == "BarrettsLike.Sig17+" | expr.matrix$status == "Sig17-",]
ERBB2_expr <- subset(samples, select = c("status", "ERBB2"))
ERBB2_expr$ERBB2 <- as.numeric(ERBB2_expr$ERBB2)

my_comparisons <- list(c("NaiveLike.Sig17+", "Sig17-"),
                       c("BarrettsLike.Sig17+", "Sig17-"))

my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3])

ERBB2_expr$status <- factor(ERBB2_expr$status, levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+","Sig17-"))

medians <- ERBB2_expr %>%
  group_by(status) %>%
  summarize(median_value = median(ERBB2, na.rm = TRUE))

pdf(file = "ERBB2_comp.pdf")
ggviolin(ERBB2_expr, x = "status", y = "ERBB2",
         fill = "status", palette = my_colors,
         add = "boxplot", add.params = list(width = 0.1, color = "black", fill = "white"))+
  xlab("")+ ylab("Expression Value") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) + geom_pwc(method = "dunn_test", label = "p.adj.format",
                                                   bracket.nudge.y = 0.5, p.adjust.method = "none",
                                                   step.increase = 0.3)
dev.off()

