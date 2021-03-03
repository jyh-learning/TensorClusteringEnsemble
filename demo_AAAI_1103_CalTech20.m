%% 5
clear
clc
addpath(genpath('.\'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataName = 'Caltech20';

members = [];
gt = [];
load(['bc_pool_',dataName,'.mat'],'members','gt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% only for CalTech %%%%%%%%%%%%
for i=1:length(gt)
    if gt(i)>18
        gt(i)=gt(i)-2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%



[N, poolSize] = size(members);
trueK = numel(unique(gt));
%%
% Number of base clusterings
M = 10;
% Repeat the proposed algorithm 
cntTimes = 20; 


% For each run, M base clusterings will be randomly drawn from the pool.
% Each row in bcIdx corresponds to an ensemble of M base clusterings.
bcIdx = zeros(cntTimes, M);
for i = 1:cntTimes
    tmp = randperm(poolSize);
    bcIdx(i,:) = tmp(1:M);
end
parfor runIdx = 1:cntTimes
%     disp('**************************************************************');
%     disp(['Run ', num2str(runIdx),':']);
%     disp('**************************************************************');
    
    %% Construct the ensemble of M base clusterings
    % baseCls is an N x M matrix, each row being a base clustering.
    baseCls{runIdx} = members(:,bcIdx(runIdx,:));
    MCA_ML{runIdx}=compute_MCA_jyh(baseCls{runIdx}) % MCA_ML denotes the coherent matrix
    CA{runIdx}=compute_CA_jyh(baseCls{runIdx});
    [A{runIdx},E{runIdx},B{runIdx}]=TensorEnsemble(MCA_ML{runIdx},CA{runIdx},.002);
end



%% For the SC
addpath(genpath('.\auxiliaryfunction'))
gnd=gt;
nClass=max(gnd);


parfor i = 1:cntTimes
[ACC_SC(i),NMI_SC(i)]=baseline_clustering_method_csc(A{i}(:,:,2)+A{i}(:,:,2)',nClass,gnd);
end

disp('Ours SC: mean ACC2 and NMI2')
disp([mean(ACC_SC),mean(NMI_SC)])
%% for EA

% for i=1:cntTimes
%    
%     resultsEA2{i} = runLWEA(A{i}(:,:,2)+A{i}(:,:,2), nClass);
%     [ACC_EA(i),NMI_EA(i),res_EA{i}]=ClusterRes2metric(resultsEA2{i},gnd);
%     
% end
% 
% disp('Ours EA: mean ACC and NMI')
% disp([mean(ACC_EA),mean(NMI_EA)])