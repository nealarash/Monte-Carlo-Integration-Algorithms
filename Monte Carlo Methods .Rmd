---
title: "Pstat 134 Final Homework 4"
author: "Neala Rashidfarrukhi"
date: "7/21/2022"
output:
  pdf_document: default
  html_document: default
  word_document: default
---

Question 1

* Gaussian Copula

Here we have  two random variables $X1-Exp(1)$ and $X2-Exp(3)$. The values 2 and 3 are rate parameters. The correlation coefficient of them is $P(X1, X2) = 0.7$.

(a) 

The question asks us to generate 100 pairs of samples under each case based on our algorithm. First we will be defining the mean of our dimension vector. Then we need to create the covariance matrix. In this case the correlation coefficient of $P(X1, X2) = 0.7$. We compute the diagonal since our covariance matrix has 0.7 in off diagonal terms since we want $Cov(Xi,Xj)=0.7$ for all i,j.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
mu = rep(0,100)
covmat <- matrix(0.7,nrow = 100,ncol =100)
diag(covmat) <- 1
```

Now we will be setting up our 100 dimensional Gaussian vector. In order to do this we must write a function for n dimensional Gaussian vector generation. In  this function n_dim_gaussvec we have three variables: dimen, mean_vec and cov_matrix. We use the t() in order to calculate the transpose of our covariance matrix. Within this t function we have the function chol() which is used for Choleski decomposition. This function computes the Choleski factorization of our covariance matrix. Then we are creating a vector which contants our transposed Choleski factorization covariance matrix and we will multiply this matrix by another matrix which contains random normally distributed variables of our dimen plus our mean vector.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
#  function for n dimensional gaussian vector generation
n_dim_gaussvec <- function(dimen, mean_vec,cov_matrix){
A<-t(chol(cov_matrix))
vec <- A %*% rnorm(dimen)+mean_vec
return(vec)
}
#  our 100 dimensional gaussian vector
gauss_vec <- n_dim_gaussvec(100,mu,covmat)
```

We are now computing the CDF of each sample. We are working with the assumption that each component is normal mean 0 and variance 1 due to the mean vector and diagonal of the covariance matrix. We are computing the CDF of each sample by using the pnorm function which gives the distribution function of our gauss_vec which we created above.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
unif_samples<- pnorm(gauss_vec)
```

Lastly I am going to be  using qexp in order to obtain correlated 100 samples of $Exp(2)$ distribution and $Exp(3)$ distribution. We are using the same code for both $X1$ and $X2$. We are using the qexp function to return the corresponding values of the quantile function. Within this function we are inputting our CDF of each sample with the rate of parameters which in this case it is $X1-Exp(2)$ and $X2-Exp(3)$

```{r, cache=TRUE, warning=FALSE, message=FALSE}
X1<- qexp(unif_samples, rate = 2)
X2<- qexp(unif_samples, rate = 3)
head(X1)
head(X2)
```

(b)

Here we will be making plots of $X2$ against $X1$. We are making a simple scatter plot using the ggplot function. As we can see there is a perfect linear relationship between the two variables.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
ggplot() +geom_point(aes(x=X1, y=X2))
```

Furthermore let's look at histograms of the $X1$ and $X2$ and see them side by side. We will be using the qplot function in r and graphing two histograms side by side. From this we can see that the two have very similar histograms. They both follow very similar trends and have a similar distribution shape. 

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot1 <- qplot(X1, color='red', main='X1~Exp(2)')
plot2 <- qplot(X2, main='X2~Exp(3)')
grid.arrange(plot1, plot2, ncol=2)
```

Here we are looking at the densities of $X1-Exp(2)$ and $X2-Exp(3)$. From this we can see $X1$ represented by the Blue density line and $X2$ represented by the Red density line. We used the ggplot function in order to create this density graph. From this we know that the density curve shoes the probability. As we can see these are not normal density curves. The density curve for $X1$ is slightly symmetrical with a slight right skew. Whereas the density curve for $X2$ is extremely right skewed. Now, let's apply a power transformation to this data.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
ggplot() + geom_density(aes(X1), color="Blue") +  geom_density(aes(X2), color="Red") 
```

Now we will be applying a power transformation to symmetrize the data.n In order to do this we will be using the boxCox function. 
The boxcox() function in MASS package computes and optionally plots profile log-likelihoods for the parameter of the Box-Cox power transformation. We must remember that this function only works for a non-ngeative response variable. In this case $X2>0$. 

First Let's verify that $X2$ >0. When using the View() I can see that there are no negative values within our response variable so we can proceed with the boxcox power transformation. Here we are taking the linear model between our response (X2) and our predictor (X1). From this we get our log likelihood. We can see that this transformation created a graph that follows the exponential distribution. Let's see how this transformed our data in other ways.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
# View(X2)
bcTrans<- boxcox(lm(X2~X1), lambda = seq(-1,1,by = .1))
```

As we can see our power transformation has transformed both of our variables but not seemingly in the way we want. The power transformation did not symmetrize the data, it looked much better before the power transformation.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot1 <- qplot(bcTrans$x, color='red', main='X1~Exp(2)')
plot2 <- qplot(bcTrans$y, main='X2~Exp(3)')
grid.arrange(plot1, plot2, ncol=2)
```

Now we will look at a simple scatter plot using the ggplot function of our power transformation. We can see that it is the same as our log-liklihood from our power transformation. 

```{r, cache=TRUE, warning=FALSE, message=FALSE}
options(scipen=0)
ggplot() + geom_point(aes(x=bcTrans$x, y=bcTrans$y), col="Blue")
```

Lastly, let's look at the linear model of the two, to see the relationship between the predictor and response before and after the power transformation. Here we are looking at the the linear_model before the transformation. We can see that our response variable has a p-value of $2e^{-16}$. Thus the relationship between the two variables has statistical significance. Furthermore, we have an $R^2$ value of 1 which is a perfect value. This means that our predictors fit our observed values. In this case, the power transformation is not necessary.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
linear_model = lm(X2~X1)
summary(linear_model)
```

Here we are looking at the predictor and response variables after the power transformation. When running this code we still have a very good p-value of our predictor value of 4.3e-07. Although, our $R^2$ value is much smaller in this case at 0.231. Thus our predictors do not fully fit our observed values. 

```{r, cache=TRUE, warning=FALSE, message=FALSE}
linear_model = lm(bcTrans$y~bcTrans$x)
summary(linear_model)
```

Question 2

* Monte Carlo Stimulation

We have two quarters of circles that are overlapped within a square. We will be using Monte Carlo stimulation in order approximate the area under the curve. We know area of this quarter circle is $\frac{\pi}{4}$. We can estimate pi by throwing lot of darts on this square with quarter circle inside. 

$$\operatorname{Pr}[\text {the dot is in the quarter circle}]=\frac{\text{area(quarter circle})}{\text {area(square) }}=\frac{\pi}{4} \text {.}$$

Define indicator variables (Bernoulli process)
$I_{i}=I\{i$-th dot is in the quarter circle $\}$, so that $I_{i} \sim_{i i d}$ Bernoulli $\left(\frac{\pi}{4}\right)$, then again, by LLN,
$$
\frac{1}{n} \sum_{i=1}^{n} I_{i} \stackrel{\mathbb{P}}{\rightarrow} \mathbb{E} I_{1}=\frac{\pi}{4}
$$

Generate a pair of (x,y) randomly by setting (x,y) = ($u_1,u_2$) where $U_1,U_2\sim Unif[0,1]$. 

**Check if $\sqrt{x_1^2+x_2^2} \leq 1$** 

The Bernoulli distribution has probability as its expectation. We are stimulating two values (u1 & u2) with our runif function. We have stimulating 1000 values for each. Since we want to have the circle within the unit circle we are using the formula from above with our stimulated random variables from the uniform distribution. Now we are creating an indicator function called est which is using an ifelse statement where our for our logical expression r<=1, then u1 is 1 and u2 is 0. From this we will take this calculated amount and divide it by 1000 which is our number of samples.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
u1 <- runif(1000)
u2 <- runif(1000)
r <-sqrt(u1**2+u2**2)
# indicator function.
est <- ifelse(r<=1,1,0)
# sum of indicator functions divided by number of samples
print(sum(est)/1000)
```

After running our algorithm we get a value. We can notice since it is random generation every time we run it we get a different value which is less than 1. Here we are getting an approximated value of 0.788. Now let's see how similar this is to our calculated quarter of a circle. 

```{r, cache=TRUE, warning=FALSE, message=FALSE}
print(pi/4)
```

When printing pi/4 we get .7854. So here we are getting a .0026 difference between our stimulated value and our calculated value. There is about a .33 % difference between the two values. Thus we can conclude that our monte carlo stimulation to approximate the shaded region was successful.

Question 3

* M-H & Gibbs Sampler

(a)

We know that that the Metropolis Hastings (M-H) algorithm is used for producing samples from distributions that may be difficult to sample from otherwise. We will stimulate a markov chain where our stationary distribution is our target distribution ($\pi$). Since we want to sample from $q(X)=C p(x)$ but we only have $p(x)$ and $C$is unknown we must sample from $p(x)$ since it is impossible to sample from $q(x)$. 

Here we are looking a t a tri-variate random vector with a given joint density. Our joint density function is defined as, 

$f(x,y,z)=C e^{-(x+y+z+xy+yz+xz)}, x,y,z>0$


In order to estimate ${\mathbb{E}}XYZ$ we will create a target distribution. Here we are using an ifelse statement that states if x<0, and our corresponding x is 0 and y is $e^x$. Now we are creating a function called MHsim that takes in 4 variables: target density, number of iterations, starting point on our MCMC and then the standard deviation of our Gaussian transition kernel. Now, it states in the question to set 500 burn-in period. The Burn-in period is when we throw away some iterations at the beginning of our Metropolis-Hastings Algorithm. Thus we are having our starting point to be 500, but we will have 10000 iterations. So in our for loop we are initializing i to go from 2 to the number of iterations we set when using our MHsim function. We have our current_x value which takes i-1, then our proposed_x which is a random simulation of the normal distribution where our probability is 1, our mean is current_x that we defined in this for loop and our standard deviation is the proposal_sd which is a input variable in our function. Now using all our created values we are creating a variable r_n where the target density of our proposed_x is divided by the target density of our current_x which is our mean. Now in our function we have an if statement that is setting the condition if r_n is higher than the $unif[0,1]$ then we will move to a different state (proposed_x) and if it is lower than the $unif[0,1]$ we will stay with our current_x. 


```{r, cache=TRUE, warning=FALSE, message=FALSE}
# target density
set.seed(1)
target_dist = function(x){
return(ifelse(x<0,0,exp(-x)))
}
MHsim<- function(target_density,niter,start_pt, proposal_sd){
x = rep(0,niter)
x[500] = start_pt
for(i in 2:niter){
current_x = x[i-1]
proposed_x = rnorm(1,mean=current_x,sd=proposal_sd)
r_n = target_density(proposed_x)/target_density(current_x)
if(runif(1)<r_n){
x[i] = proposed_x
} else {
# else, stay
x[i] = current_x
}
}
return(x)
}
```

Now we are defining the variables without our Metropolis-Hastings Algorithm algorithm where our target density is the target_dist function we created with our ifelse statement, our number of iterations is 10,000, our starting point is 500 due to our burn-in period and our proposal standard deviation increases by a multiple of 10 for each variable.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
x<- MHsim(target_dist,10000,500,1)
y<- MHsim(target_dist,10000,500,10)
z<- MHsim(target_dist,10000,500,100)
# r_n is acceptance probability
```

* Now let's visualize the data

We are going to make 3 plots on two pages. The first three plots are line graphs of each variable where we have specified in the main how our standard deviation is increasing. Now for our other three graphs we are creating histograms with a density line that is shown in red.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
par(mfrow = c(2,3))
plot(x,main="values of x visited, sd = 1",type ="l")
plot(y,main="values of x visited, sd = 10",type ="l")
plot(z,main="values of x visited, sd = 100",type ="l")
hist(x,probability = TRUE, main="Histogram of values of x, sd = 1 ")
xx = seq(0,10,length=100)
lines(xx,target_dist(xx),col="red")
hist(y,probability = TRUE, main="Histogram of values of x , sd = 10")
lines(xx,target_dist(xx),col="red")
hist(z,probability = TRUE, main="Histogram of values of x, sd = 100")
lines(xx,target_dist(xx),col="red")
```

From this we can see that as the variance of the transition kernel increases our sampling becomes less accurate.

* Here are a few other ways to look at our sampling. 

These are similar to the poiunts above we just have our graphs enlarged so they are easier for the viewer to see.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
par(mfcol=c(3,1)) #rather odd command tells R to put 3 graphs on a single page
maxz=max(c(x,y,z))
hist(x,breaks=seq(0,maxz,length=10))
hist(y,breaks=seq(0,maxz,length=10))
hist(z,breaks=seq(0,maxz,length=10))
```

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot(x,type="l")
lines(y,col=2)
lines(z,col=3)
legend(x="topright", legend=c("x","y","z"), col=c(1,2,3), lwd=2)
```

(b)

* Gibbs Sampling

Gibbs sampling is another form of MCMC methods and it is a subgroup of Metropolis-Hastings algorithm. We use Gibbs sampling when the joint distribution is not know explicitly or is difficult to sample from directly, just like the Metropolis-Hastings method. Since we know the conditional distribution we will sample from this. This sampling covers high dimensional sampling.

First we must obtain the full set of conditional distributions before sampling. Here our joint density is,

$f(x,y,z)=C e^{-(x+y+z+xy+yz+xz)}, x,y,z>0$

When calculating the marginal distributions of f(x),f(y),f(z) I used to the integral constraints $-\infty$ to $\infty$. Let's define our conditional distributions as,

$f(x,y,z) = e^{-(x+y+z+xy+yz+xz)}$
    $e^{-(x+y(1+x)+z(1+x+y))}$
    $e^{(-x)}e^{(-y(1+x))}e^{-z(1+x+y)}$
    
$f(x)=Exp(1)$
$f(y|x)=Exp(1+x)$
$f(z|x,y)=Exp(1+x+y)$

So now that we have our conditional distributions let's create our algorithm. First we are resetting everything so our sampling method is completely random. Next we are creating our function using all of our conditional distributions we found above. From this we are creating a matrix of our observations where we are initializing the matrix to have 1e6 rows and 3 columns as our dimensions. Now we are initializing a variable exp where in our for loop we are creating an observation from our function above and using the prod function to return all the values that run trough the loop

```{r, cache=TRUE, warning=FALSE, message=FALSE}
set.seed(1)
n_iteration <- 10000

gibbs_sim <- function(){
x <- rexp(1)
y <- rexp(1,1+x)
z <- rexp(1,1+x+y)
return(c(x,y,z))
}

obs <- matrix(0,n_iteration,3)
exp <- 0
for(i in 1:n_iteration){
obs[i,] <- gibbs_sim()
exp<- exp + prod(obs[i,])
}

exp/n_iteration
```

Here we have our estimate of ${\mathbb{E}}XYZ$ using gibbs sampling method.

Question 4

Here we will be using SIR (Simple Importance Resampling) to generate a permutation $1,2,...100$ whose distribution is approximately that of $X1,X2,...X100$. Here we are resetting our seed so our sampling is compltleey random. We are creating three variables that we are going to use throughout our algorithm for the permutation. We know that we need to arrange our values in an order to represent $X1,X2,...X100$.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
set.seed(1)
start<- 1:100
weight<- 1:100
set<-list()
```

Here we are creating a repeat loop where we initialize a vector. We have a while loop which takes the sum of our vector times our weight variable. When the sum times weight is less than 285000 we will sample from our start variable, 100 times with no replacement. These values represent the amount we want from our permutation. Now we are taking our set variable that we defined above to be a blank list and outside the while loop we are having the repeat loop increment through the length of each set. Then we have the condition to stop the repeat loop with our if statement.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
repeat{
  vector=0
  while(sum(vector*weight)<285000){
    vector<-sample(start,100,replace=FALSE)
  }
  set[[length(set)+1]]<-vector
  if(length(set)==1000)
    {break}
}
```

Now that we have our set list from the repeat loop we are going to use a for loop to get all of our samples into a permutation that follows $X1,X2,...X100$. So we are creating a variable collect, which uses the rep function which greats a generic function of replication, with 1000 samples. Now we have a for loop where we initialize i to run through 1 to 100. We are creating a sample variable index which samples from the values 1-1000 for only 1 value. Now we are reassigning our collect variable to take the values from our collect variable we created before and combining it with our set list of our index sample by rows.


```{r, cache=TRUE, warning=FALSE, message=FALSE}
collect<- rep(0,1000)

for(i in 1:100){
  index<- sample(1:1000,1)
  collect<-rbind(collect,set[[index]])
}
```

Now that we have created our permutation, let's look at the statistics of it. Here we have very large variances of the two variables since our data is spread out due to the sampling method.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
X1<-collect[-1,1]
X100<- collect[-1,100]
summary(X1)
summary(X100)
var(X1)
var(X100)
```

* Let's Vizualize

Here we are looking at a simple histogram of the two variables. As we can see it does not follow any specific kind of distribution. Let's see what happens if we transform the data using log and sqrt transformations.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
par(mfcol=c(2,1))
hist(X1)
hist(X100)
```
As you can see when we apply a logarithmic transformation the data is even more left skewed than it was previously. In our graphs above we can see that X1 was slgihtly right skewed but now the data has become completely left skewed. The X100 data remains unchanged just with a slight more extreme left skewness.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot1 <- qplot(log(X1), color='red', main='X1')
plot2 <- qplot(log(X100), main='X100')
grid.arrange(plot1, plot2, ncol=2)
```
Here we are trying a sqrt transformation. Once applying this transformation we can see that that our X1 data follows no distribution. If we want we could argue that the distirbution is unimodal. Furthermore as we can see on the right our X100 graph remains mostly unchanged, with more values distirbuted on the right of the graph.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot1 <- qplot(sqrt(X1), color='red', main='X1')
plot2 <- qplot(sqrt(X100), main='X100')
grid.arrange(plot1, plot2, ncol=2)
```
Lastly, to vizualize the data even further we are looking at a time series plot.

```{r, cache=TRUE, warning=FALSE, message=FALSE}
plot(X1,type="l",col=2, main="Time Series Plot of X1 & X100")
lines(X100,col=3)
legend(x="topright",legend=c("X1","X100"), col=c(2,3), lwd=2)
```


