# Study-and-Implemantation-for-Convex-and-Non-Convex-Optimization-Problems

Hi, I'm a studant at UNICAMP (Universidade Estadual de Campinas), and I
am learning about optimization techiniques, so I pretend to implement 
some of the techiniques I will learn here, with academic purpose only. The implementations here 
aren't mine but from lectures and books I find interesting. I will not 
try to do the fastest algorithms but make them easy to understand and learn.
I'm trying to comment them a lot too, I hope this helps someone :).


## Line Search

I'm just starting to implement one simple backtracking algorithm and the
zoom algorithm of the book Numerical Optimization from Jorge Nocedal and
Stephen J. Wright (second edition) on pages 60 and 61. I have to improve
the way the step direction is choose and test the algorithm better because
I was having some numerical issues.

In the Backtranking.m one can find an simple implementation of backtracking 
with steepest gradient descent

There is an implementation of the zoom algorithm (satisfyes wolf conditions
on ever iteraction) in the WolfeConditions.m file (there is a backtragking 
too, just because I was comparing the results).

In the BacktrackingWithModifiedHessian.m I compared the results with the
solution of steepest gradient descent vs a newton with modified hessian,
in this case I just sum a diagonal matrix with the value of the most negative
eigenvalue of the hessian.


