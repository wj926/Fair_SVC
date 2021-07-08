function [inv_C] = get_inv_C(expX, input);
% Compute inverse of a covariance matrix C

[n, D] = size(input); % number of training cases and dimension of input space.

% first, we write out the covariance matrix C for the training inputs ...
C = zeros(n,n);  % create and zero space
for d = 1:D      % non-linear contribution
  C = C + expX(d)*(repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2;
end
C = expX(D+1)*exp(-0.5*C) + expX(D+2)*eye(n);

inv_C = eye(n)/C;
