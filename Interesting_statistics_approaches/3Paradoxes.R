############################
#Three Bewiching Paradoxes
#Date: 02/02/2020
#Author: Zihuan Qiao
############################

# Paradox 1: two boys
# simulate 10000 two children samples
set.seed(123)
child1 <- runif(10000)
child2 <- runif(10000)
child1 <- ifelse(child1>.5, 1, 0) # 0: boy, 1: girl
child2 <- ifelse(child2>.5, 1, 0) 

# Example 4.26
ind_1st_boy <- child1==0 # first child is boy
ind_2nd_boy <- child1+child2==0 # second child is also boy
sum(ind_2nd_boy)/sum(ind_1st_boy)

# Example 4.27
ind_boy <- child1+child2!=2 # at least one boy
ind_2boys <-child1+child2==0 # both children are boys
sum(ind_2boys)/sum(ind_boy)

# Paradox 2: Switch envelope
# Example 4.28
# simulate 10000 Ali & Baba envelope cases
set.seed(123)
r <- runif(10000)
Ali <- ifelse(r>.5, 10, 5) 
Baba <- ifelse(r<.5, 10, 5) 
mean(Ali)
mean(Baba)

# Paradox 3: Switch envelope strategy
switch_envelope <- function(rate){
  U <- rexp(10000, rate = rate)
  Ali_switch <- ifelse(Ali<U, Baba, Ali)
  Baba_switch <- ifelse(Ali<U, Ali, Baba)
  mean(Ali_switch>Baba_switch)
}
exp_rate <- seq(.05, 1, .05)
Ali_gain <- sapply(exp_rate, switch_envelope)
plot(x=exp_rate, y=Ali_gain, xlab = "exponential rate", ylab = "prob of Ali>Baba", type = "l")
abline(h=mean(Ali>Baba), col=2)

