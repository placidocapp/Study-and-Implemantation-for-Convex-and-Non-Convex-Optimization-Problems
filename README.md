# Study-and-Implemantation-for-Convex-and-Non-Convex-Optimization-Problems

Hi, I'm a studant at UNICAMP (Universidade Estadual de Campinas), and I
am learning about optimization techiniques, so I pretend to implement 
some of the techiniques I will learn here, with academic purpose only. The implementations here 
aren't mine but from lectures and books I find interesting. I will not 
try to do the fastest algorithms but make them easy to understand and learn.
I'm trying to comment them a lot too, I hope this helps someone :).


## Line Search

Here you will find the newton method and steepest gradient in one file and
the quasi newton methods in another:

#### Newton and Steepest

* backtrack.m file contains the famous backtraking algorithm
* lineSearch.m is a function to choose the step size, it calls the zoom.m function that's a function
describled in Numercial Optimization book that satisfies the wolfe conditions
* mcfac.m contains a modified cholesky factorization to correct the hessian in case it is not positive definite
* ModifiedNewtonAndSteepest.m contains the newton and steepest methods with the linesearch step and the cholesky modifications
* SteepestDescentBacktracking.m here is the same as the last function exept that the step is choose by the backtracking algorithm
in general backtracking is faster in each step but it's performance is worse than lineSearch method comparing the number of
steps.


