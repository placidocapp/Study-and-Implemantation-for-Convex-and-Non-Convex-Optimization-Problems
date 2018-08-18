# Study-and-Implemantation-for-Convex-and-Non-Convex-Optimization-Problems

Hi, I'm a student at UNICAMP (Universidade Estadual de Campinas), and I
am learning about optimization techniques, so I pretend to implement 
some of the techniques I will learn here, with academic purpose only. The implementations here 
aren't mine but from lectures and books I find interesting. I will not 
try to do the fastest algorithms but make them easy to understand and learn.
I'm trying to comment them a lot too, I hope this helps someone :).


## Line Search

Here you will find the newton method and steepest gradient in one file and
the quasi newton methods in another:

#### Newton and Steepest

* backtrack.m file contains the famous backtraking algorithm
* lineSearch.m is a function to choose the step size, it calls the zoom.m function that's a function
describled in Numerical Optimization book that satisfies the wolfe conditions
* mcfac.m contains a modified cholesky factorization to correct the hessian in case it is not positive definite
* ModifiedNewtonAndSteepest.m contains the newton and steepest methods with the line search step and the cholesky modifications
* SteepestDescentBacktracking.m here is the same as the last function except that the step is choose by the backtracking algorithm
in general backtracking is faster in each step but it's performance is worse than lineSearch method comparing the number of
steps.
*zoom.m is an algorithm proposed in the book of Numerical Optimization to find an step satisfying the Wolfe conditions.

#### Quasi-Newton

*QuasiNewtonComparison.m and QuasiNewtonMethods.m are quite similar, the only difference is that in the first one all the 
algorithms solve the problem and in the second one you can choose your algorithm and play. The objective is to make it 
easy to compare the methods but if one want to solve a problem the second is probably more appropriate.
*SR1.m, the book of numerical optimizations suggests an trust region method to SR1 so I programmed it
*bisec.m is a method to find an step length that satisfies the wolfe conditions, However I tested it and the zoom algorithm 
performs much better

## Trust Region

*In the TrustRegion.m we have al the methods so you can choose one of them and run
*TrustRegion.m just make it easy to compare the algorithms

## Gradient

I still have to check if everything here is correct, I had a little doubt in some methods. Sometimes the step of the 
gradient becomes very small and the algorithms did not converge yet so maybe I am doing something wrong.

#### Logistic Regression Just for Train

*Here I made a logistic regression to given a number, find out if it is a pre-selected number or not, I made it just
to study and learn the logistic regression before implementing the neural network itself.

#### Others

*NeuralNetwork.m Here I implement all the gradients and run them. Note that the algorithm gets a lot slower because I have
to calculate the function cost for the training set every iteration to show a graph of the evolution of each method.
I am practice no one would calculate it every time, that slows the algorithm a lot.
*g.m and dg.m are the activation functions and its derivatives, respectively.
*fgrad.m do the backpropagation to calculate the derivatives of the weights on the neural network.
*fcost.m runs every iteration to calculate the value of objective function. (for a batch of the training set)

## Linear Optimization

*Working on that ...
