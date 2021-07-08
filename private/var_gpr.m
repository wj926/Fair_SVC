function var = var_gpr(test, input, inv_C, hyperpara)
% Variance in Gaussian Process Regression used as Support Funtion of Clustering
%
% The variance function of a predictive distribution of GPR
% sigma^2(x) = kappa - k'C^(-1)k
%==========================================================================
% Implemented by H.C. Kim Jan. 16, 2006
% Modified by Sujee Lee at September 10, 2014.
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

[n,D]=size(input);
[nn,D]=size(test);

expX=hyperpara;
% X=log(expX);

a = zeros(n, nn);  % create and zero space

for d = 1:D
  a = a + expX(d)*(repmat(input(:,d),1,nn)-repmat(test(:,d)',n,1)).^2;
end
a = expX(D+1)*exp(-0.5*a);
b = expX(D+1);

s_a_inv_C_a = sum(a.*(inv_C*a),1);
var = b - s_a_inv_C_a';
%va_deri=0;