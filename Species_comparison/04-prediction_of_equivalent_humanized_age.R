library("tidyverse")
library("SummarizedExperiment")
library("glmnet")
library("lspline")
library("scales")

species_comparison <- readRDS("Results/species_comparison.rds")
human_comparison <- species_comparison[, species_comparison$species == "human"]
x <- t(assay(human_comparison))
y <- human_comparison$log2_age_days

#compare between different values of alpha
foldid <- sample(1:10, size = length(y), replace = TRUE)
cv1  <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), foldid = foldid, alpha = 1)
cv.5 <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), foldid = foldid, alpha = 0.5)
cv.25 <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), foldid = foldid, alpha = 0.25)
cv.1 <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), foldid = foldid, alpha = 0.1)
cv0  <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), foldid = foldid, alpha = 0)

par(mfrow = c(2,3))
plot(cv1); plot(cv.5); plot(cv.25); plot(cv.1); plot(cv0)
plot(log(cv1$lambda)   , cv1$cvm , pch = 19, col = "red",
     xlab = "log(Lambda)", ylab = cv1$name)
points(log(cv.5$lambda) , cv.5$cvm , pch = 19, col = "orange")
points(log(cv.25$lambda) , cv.25$cvm , pch = 19, col = "grey")
points(log(cv.1$lambda), cv.1$cvm, pch = 19, col = "blue")
points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "green")
dev.off()
#choose alpha = 0 and do ridge regression
#use cross validation to choose the optimal lambda
set.seed(1)
cvfit <- cv.glmnet(x, y, lambda=10^(seq(-4,7,0.05)), alpha = 0)
plot(cvfit)
print(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#choose lambda = cvfit$lambda.min
set.seed(1)
fit <- glmnet(x, y, lambda = cvfit$lambda.min, alpha = 0)
print(fit)
coef(fit)

#predict humanized ages for mouse and macaque data
macaque_comparison <- species_comparison[, species_comparison$species == "macaque"]
macaque_x <- t(assay(macaque_comparison))
macaque_y <- predict(fit, newx = macaque_x, type = "response")
colnames(macaque_y) <- "humanized_age"
macaque_comparison$humanized_log2_age_days <- as.vector(macaque_y)
macaque_comparison$humanized_age <- 2^(macaque_comparison$humanized_log2_age_days)

mouse_comparison <- species_comparison[, species_comparison$species == "mouse"]
mouse_x <- t(assay(mouse_comparison))
mouse_y<- predict(fit, newx = mouse_x, type = "response")
mouse_comparison$humanized_log2_age_days <- as.vector(mouse_y)
mouse_comparison$humanized_age <- 2^(mouse_comparison$humanized_log2_age_days)

human_comparison <- species_comparison[, species_comparison$species == "human"]
human_x <- t(assay(human_comparison))
human_y<- predict(fit, newx = human_x, type = "response")
human_comparison$humanized_log2_age_days <- as.vector(human_y)
human_comparison$humanized_age <- 2^(human_comparison$humanized_log2_age_days)

humanized_data <- as.data.frame(rbind(cbind(human_comparison$label, human_comparison$humanized_log2_age_days, human_comparison$humanized_age),
                                 cbind(macaque_comparison$label, macaque_comparison$humanized_log2_age_days, macaque_comparison$humanized_age),
                                 cbind(mouse_comparison$label, mouse_comparison$humanized_log2_age_days, mouse_comparison$humanized_age)
                                 ))

colnames(humanized_data) <- c("label", "humanized_log2_age_days", "humanized_age")
humanized_data$humanized_log2_age_days <- as.numeric(humanized_data$humanized_log2_age_days)
humanized_data$humanized_age <- as.numeric(humanized_data$humanized_age)

#add humanized data to species_comparison
coldata <- as.data.frame(colData(species_comparison))
coldata2 <- coldata %>% left_join(humanized_data)

species_comparison$humanized_age <- coldata2$humanized_age
species_comparison$humanized_log2_age_days <- coldata2$humanized_log2_age_days

saveRDS(species_comparison, "Results/species_comparison_humanized_age.rds")


#plot age transformation
p <- ggplot(data = coldata2, mapping = aes(x = log2_age_days, y = humanized_log2_age_days, color = species)) +
  geom_point(size = 2.5, alpha = 0.8, show.legend = T) +
  coord_fixed(ratio = 1) +
  scale_y_continuous(limits = c(4, 14.5), breaks = c(6.8073549221,8.0552824355,9.30149619498,10.7532167492,12.885315061, 14.1764848473),
                     labels = c("GW18", "Birth", "Year01", "Year04", "Year20", "Year50")) +
  scale_x_continuous(limits = c(4, 14.5),breaks = c(4.24792751344,6,7.37503943135,8.0552824355,10,12,14),
                     labels = c(19, 64, 166, 266,1024,4096,16384),
                     sec.axis = dup_axis(name = NULL, breaks = c(4.24792751344,7.37503943135,8.0552824355),
                                         labels = c("Mouse birth", "Macaque birth", "Human birth"))) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "gray80", linetype = 3),
        panel.grid.major.x = element_line(colour = "gray80", linetype = 3),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black")
  ) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.text.y = element_text(face = c("plain", "bold", "plain", "plain", "plain", "plain"),hjust = c(0.5,0.5,0.5,0.5,0.5,0.5)),
        axis.text.x.top = element_text(angle = 45, hjust = 0),
        axis.title = element_text(size=14)) +
  theme(legend.key.size = unit(0.2, 'inch'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position="top") +
  geom_vline(xintercept = c(4.24792751344,7.37503943135,8.0552824355), linetype = 3, color = hue_pal()(3)[3:1]) +
  geom_hline(yintercept = 8.0552824355, linetype = 3) +
  # geom_abline(slope = 1, intercept = 0, colour = "gray60", linetype = 3) +
  labs(y = "Predicted equivalent human PSD age",
       x = "Real post-conceptional age (days)")

human_slope <- coldata2 %>% filter(species == "human")
human_slope_coef <- lm(humanized_log2_age_days ~ log2_age_days, data = human_slope)$coefficients[2]
mouse_slope <- coldata2 %>% filter(species == "mouse")
mouse_slope_coef <- lm(humanized_log2_age_days ~ log2_age_days, data = mouse_slope)$coefficients[2]
macaque_slope <- coldata2 %>% filter(species == "macaque")
macaque_slope$x <- macaque_slope$log2_age_days
macaque_slope$y <- macaque_slope$humanized_log2_age_days
m1 <- lm(y ~ lspline(x, knots = 8.366322), data=macaque_slope)
m1$coefficients

p + geom_smooth(
  data = human_slope, aes(x = log2_age_days, y = humanized_log2_age_days),  # grouping variable does the plots for us!
  method = "lm", se = FALSE, color = hue_pal()(3)[1],
  formula = y ~ x, linetype = "dashed"
) + geom_smooth(
  data = mouse_slope, aes(x = log2_age_days, y = humanized_log2_age_days),  # grouping variable does the plots for us!
  method = "lm", se = FALSE, color = hue_pal()(3)[3],
  formula = y ~ x, linetype = "dashed"
) + geom_smooth(
  data = macaque_slope, aes(x = log2_age_days, y = humanized_log2_age_days),  # grouping variable does the plots for us!
  method = "lm", se = FALSE, color = hue_pal()(3)[2],
  formula=formula(m1), linetype = "dashed"
)


ggsave("Results/species_comparison/predicted_humanized_PSD_age_across_species.pdf", width = 6, height = 6)



