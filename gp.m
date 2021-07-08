function model=gp(input,hyperparams)
% Gaussian Process Support Function for Clustering 
%
%  Gaussian process support function which is the variance function of a
%  predictive distribution of GPR : 
%   sigma^2(x) = kappa - k'C^(-1)k
%  where covariance matrix C(i,j) is a parameterized function of x(i) and
%  x(j) with hyperparameters Theta, C(i,j) = C(x(i),x(j);Theta),
%  kappa = C(x~,x~;Theta) for a new data point x~ = x(n+1),
%  k = [C(x~,x1;Theta),...,C(x~,x(n);Theta)]
%
% Synopsis:
%  model = gp(input)
%  model = gp(input,hyperparams)
%
% Description :
% It computes variance function of gaussian process regression learned from
% a training data which can be an estimate of the support of a probability
% density function. A dynamic process associated with the variance function
% can be built and applied to cluster labeling of the data points. The
% variance function estimates the support region by sigma^2 <= theta, where
% theta = max(sigma^2(x))
%
% Input:
%  input [dim x num_data] Input data.
%  hyperparams [(num_data + 2) x 1] 
%
% Output:
%  model [struct] Center of the ball in the kernel feature space:
%   .input [dim x num_data] 
%   .hyperparams [(num_data + 2) x 1]  
%   .inside_ind [1 x num_data] : all input points are required to be used for computing
%   support function value.
%   .inv_C [num_data x num_data] : inverse of a covariance matrix C for the training inputs
%   .r : support function level with which covers estimated support region 
%
%==========================================================================
% Implemented by Kyu-Hwan Jung.
% Modified by Sujee Lee at September 10, 2014.
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================
tic

disp('Step 1 : Training Support Function by GP...')

model.X = input';
model.hyperparams = hyperparams;
model.inside_ind = [1:size(model.X,1)];
model.inv_C = get_inv_C(model.hyperparams,model.X);

tmp = var_gpr(model.X, model.X, model.inv_C, model.hyperparams);

R = max(tmp);
model.r = R;

training_time_gp = toc;

disp(['Training Completed !' char(10)])