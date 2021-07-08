function [sinput, support_model] = modelsupport(input, support, options, varargin) 

default_support = 'SVDD';
default_options = struct('ker','rbf','arg',0.1,'C',0.5, 'gpparam',[100*ones(size(input,2),1); 1; 10]); 
% Parameters for SVDD :
%  ker(kernel) : 'rbf'(Gaussian) / arg : sigma in kernel / C : constant
% Parameters for GP :
%  gpparam : Hyper-parameters for GP
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

if nargin < 2
    support = default_support;
    if nargin < 3
        options = default_options;
    end
end

%% Step 1 : Training SVDD

switch support
    case 'SVDD'%==== SVDD for support function
        sinput = my_normalize(input');   
        support_model = svdd(sinput,options);
        
    case 'GP'%===== GP for support function 
        sinput = my_normalize2(input); 
        hyperparams = options.gpparam; %Hyper-parameters for GP
        support_model = gp(sinput, hyperparams);
  
end
support_model.support_type = support;