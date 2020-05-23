


###################################
sampleN = 20
nCG = 10000
trueProp = 0.05
cgProps = rep(-1, nCG)
for (i in 1:nCG) {
  cgProps[i] = mean(sample(c(0, 1), size = sampleN, replace = TRUE, prob = c(1-trueProp, trueProp)))
}
hist(cgProps, breaks = seq(0, 1, by=0.05))
