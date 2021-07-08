%==========================================================================
%   Normalizing inputs for GP
%
%   Implemented by Kyu-Hwan Jung, April 26, 2010
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

function output=my_normalize2(input)
% input [num_data x dim]
% output [dim x num_data]

[num, dim]=size(input);

max_x_tr = max(input(:)); 
min_x_tr = min(input(:));
output = (input-repmat(min_x_tr,num,dim))./repmat(max_x_tr-min_x_tr,num,dim);  

output=output';