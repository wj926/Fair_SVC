function d = kradius(X,model)

% =========================================================
%
% KRADIUS computes the squared distance between vector in kernel space 
% and the center of support. 
% 
% Implemented by Kyu-Hwan Jung
% April 26, 2010. 
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 

switch model.support_type
    case 'SVDD'

    [dim,num_data]=size(X);

    x2 = diagker( X, model.options.ker, model.options.arg);

    Ksvx = kernel( X, model.sv.X, model.options.ker, model.options.arg);

    d = x2 - 2*Ksvx*model.Alpha(:) + model.b*ones(num_data,1) ;

    case 'GP'

    for i=1:size(X,2)
       %[predict_label, accuracy] = var_gpr(X(:,i)', model.X, model.inv_C, model.hyperparams);
       predict_label = var_gpr(X(:,i)', model.X, model.inv_C, model.hyperparams);
       d(i,:)=predict_label; 
    end
        
end
