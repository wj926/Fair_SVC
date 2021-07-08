clc;
clear all;

%% Setting Input Data

load 'data/20000_adult_reduced.csv';
%load 'data/bank_sample5.csv';
%load 'data/gaussiandata.csv';

%load data/toy.mat;  
%inputwithoutsensitive = banksample(:,1:3);
inputwithoutsensitive = table(moon_not_rand(:,1:2));

%inputwithoutsensitive = table(gaussiandata(:,1:2));

finalsensitive = table(moon_not_rand(:,3));
%finalsensitive = table(gaussiandata(:,3));

finalsensitive = table2array(finalsensitive);
inputwithoutsensitive = table2array(inputwithoutsensitive);
input = cat(2,inputwithoutsensitive,finalsensitive);
%input=X';
%output=y';

%%Ring sensitive
%load 'data/ring_sensitive.mat';
%input = ring_sensitive;
%inputwithoutsensitive = input(:,1:2);
%finalsensitive = input(:,3);

%% Setting Input Data
%load 'data/adult_reduced_sample1.csv';
%load 'data/adult_sample1.csv';
%inputwithoutsensitive = table(adult_reduced_sample1);
%finalsensitive = table(adult_sample1(:,6));
%finalsensitive = table2array(finalsensitive);
%inputwithoutsensitive = table2array(inputwithoutsensitive);
%input = cat(2,inputwithoutsensitive,finalsensitive);

load 'real_data/152_job_reduced.csv';
load 'real_data/152_job_sample.csv';
inputwithoutsensitive = table(X152_job_reduced);
finalsensitive = table(X152_job_sample(:,-1));
finalsensitive = table2array(finalsensitive);
inputwithoutsensitive = table2array(inputwithoutsensitive);
input = cat(2,inputwithoutsensitive,finalsensitive);

%load 'data/adult_reduced_sample1.csv';
%load 'data/adult_sample1.csv';
%inputwithoutsensitive = table(adult_sample1(:,1:5));
%finalsensitive = table(adult_sample1(:,6));
%finalsensitive = table2array(finalsensitive);
%inputwithoutsensitive = table2array(inputwithoutsensitive);
%input = cat(2,inputwithoutsensitive,finalsensitive);
%% Two ring data
%[inputwithoutsensitive, input, finalsensitive] = makepizza(200,3,5);

%%
Result = zeros(10,5);
for KK=2:9
    method = 'T-MSC'; % method = {'CG-SC', 'S-MSC', 'T-MSC', 'F-MSC', 'V-MSC'}
    fairness = 'balanceincrease'; %fairness = {'nofair','balanceincrease','worstupdate'}
    support = 'SVDD'; % support type = {'SVDD', 'GP'}
    supportopt = struct('ker','rbf','arg',0.35,'C',0.5,... % SVDD
        'gpparam',[100*ones(size(input,2),1); 1; 10]); %GP
    options = struct('hierarchical',true,'K',KK,'epsilon',0.05,...
        'R1',0.01,'R2',0.01);
    voronoiopt = struct('samplerate', 0.2,'voronoiMethod','T-MSC',...
        'hierarchical',true,'K',4,'epsilon',0.05,... % T-SVC
        'R1',0.02,'R2',0.02,... % F-SVC
        'eta',0.1);

    [tinput, clstmodel] = msclustering(input, method, support, supportopt, options, voronoiopt,fairness);

     %% Calculating the Balance
    clusterwithsensitive = cat(1, reshape(clstmodel.cluster_labels,[1,2000]), transpose(finalsensitive));
    Nofcluster = max(clstmodel.cluster_labels);
    finalbalance = zeros(Nofcluster,1);
    for i=1:Nofcluster
        first = sum(clusterwithsensitive(2,clusterwithsensitive(1,:)==i));
        [a,total] = size(clusterwithsensitive(2,clusterwithsensitive(1,:)==i));   
        second = total - first;
        finalbalance(i) = min([first/total, second/total]);
        %fprintf('The %d th Balance is %f \n',i,finalbalance(i));
        %fprintf('Total : %d \n First : %d \n Second : %d \n',total,first,second);
    end
    %fprintf('The %d th Balance is %f \n',KK,min(finalbalance));
    
    %Results

    Result(KK,1) = KK;
    Result(KK,2) = min(finalbalance);
    %Sil = evalclusters(inputwithoutsensitive, clstmodel.cluster_labels','Silhouette');
    %Result(KK,3) = Sil.CriterionValues;
    %Calin = evalclusters(inputwithoutsensitive, clstmodel.cluster_labels','CalinskiHarabasz');
    %Result(KK,4) = Calin.CriterionValues;
    %Davies = evalclusters(inputwithoutsensitive, clstmodel.cluster_labels','DaviesBouldin');   
    %Result(KK,5) = Davies.CriterionValues;
 
    %if size(tinput,1)==2
    %    plotmsc(tinput,clstmodel);    
    %    title('T-SVC SIGMA 0.2 C 0.5');
    %end
    save = transpose(cat(1, reshape(clstmodel.cluster_labels,[1,2000]), transpose(input)));
    filename = ['./results/ours_' num2str(KK) '_adult.csv'];
    writematrix(save, filename);
end


