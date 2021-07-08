function [stan_input, clstmodel] = msclustering(input, method, support, supportopt, options, voronoiopt, fairness, varargin)
% MSCLUSTERING Multi-Basin Support-Based Clustering Method
% [stan_input, clstmodel] = svclustering(input, method, support, supportopt, options, voronoiopt)
% finds cluster labels for the input data based on support functions. 
%
% METHODS indicates the support vector clustering method, especially its
% labeling methods, which are
%     CG-SC : complete-graph based labeling method
%         (A. Ben-Hur, D. Horn, H.T. Siegelmann, and V. Vapnik, Support
%         Vector Clustering)
%     S-MSC : support function dynamics based labeling method using
%     stable-equilibrium point
%         (J. Lee and D. Lee, An Improved Cluster Labeling Method for
%         Support Vector Clustering)
%     T-MSC : transition point based adjacency decision labeling method
%         (J. Lee and D. Lee, Dynamic Characterization of Cluster
%         Structures for Robust and Inductive Support Vector Clustering) 
%     with hierarchical clustering using dynamic dissimilarity measure
%         (D. Lee and J. Lee, Dynamic Dissimilarity Measure for
%         Support-Based Clustering) 
%     F-MSC : fast support-based cluster labeling method
%         (K. Jung, D. Lee and J. Lee, Fast support-based clustering method
%         for large-scale problems) 
%     V-MSC : voronoi cell-based kernel support clustering
%         (K. Kim, Y. Son and J. Lee, Voronoi Cell-based Kernel Support
%         Clustering) 
%
% SUPPORT indicates methods constructing the support functions, which are
%     SVDD : Radius Function of Support Vector Domain Description
%         (D. Tax and R. Duin, Support Vector Domain Description)
%     GP : Gaussian process support function which is the variance function
%     of a predictive distribution of GPR 
%         (H. Kim and J. Lee, Clustering Based on Gaussian Processes)
% 
% SUPPORTOPT
%     ker : kerner type. 'rbf' type is desired
%     arg : Argument for kernel and it is different from each kernel 
%     C : Regularization constant (default, C=1)
%     hyperparam : the length scale of autocorrelation, the overall scale
%     of the function, and the amount of noise (for GP)
% 
% OPTIONS 
%     hierarchical : Use dynamic dissimilarity to do hierarchical clustering 
%     K : desired number of clusters when using the hierarchical labeling
%     epsilon : small value for moving the starting point in finding TPs
%     R1 : initial radius of balls in partition phase of the fast support-based clustering
%     R2 : intermediate ball merging parameter in the fast support-based clustering
%
% VORONOIOPT 
%     samplerate : data sampling rate for estimate support function and labeling 
%     voronoiMethod : labeling method for sample data in voronoi cell based clustering
% (Options of the labeling method(voronoiMethod) for the sample data : )
%     hierarchical : Use dynamic dissimilarity to do hierarchical
%     clustering 
%     K : desired number of clusters when using the hierarchical labeling
%     epsilon : small value for moving the starting point in finding TPs
%     R1 : radius of balls in partition phase of the fast support-based clustering
%     R2 : ball merging parameter in the fast support-based clustering
%     rho : fixed constant to compare with the ambiguity ratio for slow phase
%==========================================================================
%
% * The source code is available under the GNU LESSER GENERAL PUBLIC
% LICENSE, version 2.1. 
%==========================================================================

default_method = 'T-MSC';
default_support  = 'SVDD';
default_supportopt = struct('ker','rbf','arg',0.1,'C',0.5,... % SVDD
    'gpparam',[100*ones(size(input,2),1); 1; 10]); %GP
default_opt = struct('hierarchical',false,'K',4,'epsilon',0.05,... % T-SVC
    'R1',0,'R2',0); % F-SVC
default_vrnopt = struct('samplerate', 0.1,'voronoiMethod','F-MSC',...
    'hierarchical',false,'K',3,'epsilon',0.05,... % T-SVC
    'R1',0.01,'R2',0.01,... % F-SVC
    'rho',0.1);

if nargin < 6
    voronoiopt = default_vrnopt;
    if nargin < 5
        options = default_opt;
        if nargin < 4
            supportopt = default_supportopt;
            if nargin < 3
                support = default_support;
                if nargin < 2
                    method = default_method;
                    if nargin < 1
                        input = [];
                    end
                end
            end
        end
    end
end

switch method
    case 'CG-SC'
        [stan_input, clstmodel] = cgsc(input, support, supportopt);          
    case 'S-MSC'
        [stan_input, clstmodel] = smsc(input, support, supportopt);       
    case 'T-MSC'
        [stan_input, clstmodel] = tmsc(input, support, supportopt, options, fairness);       
    case 'F-MSC'
       [stan_input, clstmodel] = fmsc(input, support, supportopt, options);     
    case 'V-MSC'
       [stan_input, clstmodel] = vmsc(input, support, supportopt, voronoiopt);     
end
clstmodel.labeling_method = method;

  