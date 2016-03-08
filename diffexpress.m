%% Start of the mathworks workflow
load melRNAseq.mat;
melData(1:10,:)

% Took a peek at the data, want to assign groups for ablated and unablated
ablated=melFPKM(1:6,:);
unablated=melFPKM(7:12,:);

% Workflow suggests that we can make a summary table of reads that were
% assigned vs unassigned if our files came with that information.

% Chromosome mapping portion is dependent on having chromosome position
% information in the starting data set... Shoot.

% Also, need read count data to normalize and everything that follows from
% that.