# short examples of mathematical concepts

set.seed(1234)

# show that linear regression gives smaller coefficients to correlated variables
# while PCA does not
a = -5:5
b = a + rnorm(11) 
plot(a, b)
d = c(1:5, -3:-1, 2:4)
depVar = a + rnorm(11)
lm(formula= depVar ~ a + b + d)
# a and b are highly correlated
cor.test(a, b)
# remove a from model and coefficient for b increases greatly
# d increases just a little
lm(formula= depVar ~ b + d)
lm(formula= depVar ~ a + d)

# regression is like a zero sum game: there is limited variation
# to explain and the ability to explain the variation
# must be split among the independent variables
# for PCA, when you add more correlated variables, 
# there is simply more variation to be explained
veryCorMat = matrix(rep(1:5, 5), 5, 5)
prcomp(veryCorMat, center = TRUE, scale. = FALSE)
# when adding a correlated variable, the overall variation in that direction 
# increases. The loadings do decrease but
# no correlated variable is preferred over the others (they stay equal)
veryCorMat = matrix(rep(1:5, 6), nrow = 5, ncol = 6)
prcomp(veryCorMat, center = TRUE, scale. = FALSE)

veryCorMat = matrix(rep(1:5, 100), nrow = 5, ncol = 100)
prcomp(veryCorMat, center = TRUE, scale. = FALSE)$x
# the variability does not increase linearly with the number of variables added
increasingVar = function(x) dist(x = rbind(rep(1, x), rep(5, x)))
plot(1:100, sapply(X = 1:100, increasingVar))
# how is a less correlated variable affected?
veryCorMat = matrix(c(rep(1:5, 6), c(2, 1, 3, 2, 1)), nrow = 5, ncol = 7)
prcomp(veryCorMat, center = TRUE, scale. = FALSE)$rotation
veryCorMat = matrix(c(rep(1:5, 100), c(2, 1, 3, 2, 1)), nrow = 5, ncol = 101)
prcomp(veryCorMat, center = TRUE, scale. = FALSE)$rotation
veryCorMat = matrix(c(rep(1:5, 1000), c(2, 1, 3, 2, 1)), nrow = 5, ncol = 1001)
tail(prcomp(veryCorMat, center = TRUE, scale. = FALSE)$rotation)
# comparing 2nd to last row to the last row, the less correlated
# variable gets smaller ~proportionally to how much the very 
# correlated variables get smaller