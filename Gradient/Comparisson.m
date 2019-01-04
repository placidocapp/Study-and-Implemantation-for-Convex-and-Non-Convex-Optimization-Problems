%% Plot comparisson for 10000 training examples

clear all

%Load the Batch gradient
load('BatchGradient.mat')
Time_Optimizing = time;
Accuracy = accuracy;
Best_Sol = best_J;
plot(J), hold on

%Load the Stochastic gradient
load('Stochastic gradient descent.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient
load('MiniBatchGradient.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with momentum
load('MiniBatchGradientWithMomentum.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithNesterov.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithAdagrad.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

Method = {'Batch Gradient'; 'Stochastic gradient';'Mini Batch Gradient';...
    'Mini Batch Gradient With Momentum';'Mini Batch Gradient With Nesterov';...
    'Mini Batch Gradient With Modified Adagrad'};
legend(Method)
table(Method, Time_Optimizing, Accuracy, Best_Sol)

%% Plot comparisson for 100 training examples

clear all
figure()

%Load the Batch gradient
load('BatchGradient_100.mat')
Time_Optimizing = time;
Accuracy = accuracy;
Best_Sol = best_J;
plot(J), axis([0 1000 0 6]), hold on

%Load the Stochastic gradient
load('StochasticGradientDescent_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient
load('MiniBatchGradient_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with momentum
load('MiniBatchGradientWithMomentum_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithNesterov_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithAdagrad_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithAdadelta_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithRMSprop_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J)

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithAdam_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J,'g')

%Load the Mini Batch gradient with nesterov
load('MiniBatchGradientWithAdamax_100.mat')
Time_Optimizing = [Time_Optimizing; time];
Accuracy = [Accuracy; accuracy];
Best_Sol = [Best_Sol; best_J];
plot(J,'k')

Method = {'Batch Gradient'; 'Stochastic gradient';'Mini Batch Gradient';...
    'Mini Batch Gradient With Momentum';'Mini Batch Gradient With Nesterov'...
    ;'Mini Batch Gradient With Adagrad';'Mini Batch Gradient With Adadelta';...
    'Mini Batch Gradient With RMSprop'; 'Mini Batch Gradient With Adam'...
    ; 'Mini Batch Gradient With Adamax'};
legend(Method)
table(Method, Time_Optimizing, Accuracy, Best_Sol)
