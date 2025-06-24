library(ggplot2)

n <- 10000

# prior modèle 1
log_prior1 <- rnorm(n, 0, 1)
prior1 <- exp(log_prior1)

# hist(prior1, nclass = 50)
# hist(prior1[prior1 < 50], nclass = 50)



# prior modèle 5
sd.b <- runif(n, 0, 5)
mu.b <- rnorm(n, 0, .5)
log_prior5 <- NULL
for (i in 1:n) {log_prior5[i] <- rnorm(1, mu.b[i], sd.b[i])}
prior5 <- exp(log_prior5)

# hist(log_prior5, nclass = 50)
# hist(prior5, nclass = 50)
# hist(prior5[prior5 < 50], nclass = 50)

data_log <- data.frame(
  type = c(rep("modèle 1", n), rep("modèle 5",n)),
  value = c(log_prior1, log_prior5)
)

p <- data_log %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("coral", "royalblue")) +
  labs(fill="")
p

data <- data.frame(
  type = c(rep("modèle 1", n), rep("modèle 5",n)),
  value = c(prior1, prior5)
)
# 
# q <- data %>%
#   ggplot( aes(x=value, fill=type)) +
#   geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
#   scale_fill_manual(values=c("coral", "royalblue")) +
#   labs(fill="")
# q


data_petit <- data %>%
  filter(value < 10)

r <- data_petit %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity', binwidth = 0.1) +
  scale_fill_manual(values=c("coral", "royalblue")) +
  labs(fill="")
r
