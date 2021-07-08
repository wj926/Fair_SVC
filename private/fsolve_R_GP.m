function [g,H] = fsolve_R_GP(x,model)
%
%   Calculating gradient and Hessian of the trained kernel radius function
%
%==========================================================================
% Implemented by Sujee Lee at September 3, 2014.
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

input = model.X;
hparam = model.hyperparams;

[N,D]=size(input);
[nn,D]=size(x); % nn=1

inv_K = model.inv_C;


k = zeros(N, nn);  % create and zero space
for d = 1:D
    k = k + hparam(d)*(repmat(input(:,d),1,nn)-repmat(x(:,d)',N,1)).^2; % N * 1
end
k = hparam(D+1)*exp(-0.5*k);

gk = - repmat(k,1,D).*repmat(hparam(1:D)',N,1).*(repmat(x,N,1)-input); % N * D
g = - 2 * gk' * inv_K * k;
%     g= gradient(@(xx)var_gpr(xx,model.X,model.inv_C,model.hyperparams),x); %%

if nargout>1
    Hk = zeros(N,D^2); % Hessian of k. (make N * D^2 matrix)
    for j = 1:D
        Hkj = gk.*(-hparam(j)*repmat((x(:,j)-input(:,j)),1,D));  % x(:,j)-input(:,j) : scalar - vector
        Hkj(:,j) = Hkj(:,j) - hparam(j)*k;
        Hk(:,((j-1)*D+1):(j*D)) = Hkj;
    end
    H1 = Hk'*inv_K*k; % D^2 * 1
    H1 = reshape(H1,D,D); % D * D
    H2 = gk'*inv_K*gk; % D * D
    H = -2*(H1+H2);
    %         H= hessian(@(xx)var_gpr(xx,model.X,model.inv_C,model.hyperparams),x); %%
end

