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

% Normalization
% estimate pseudo-reference with geometric mean row by row
melFPKMalt=melFPKM';
pseudoRefSample = geomean(melFPKMalt,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,melFPKMalt(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

%figure;

%subplot(2,1,1)
%maboxplot(log2(melFPKMalt),'title','Raw Read Count','orientation','horizontal')
%ylabel('sample')
%xlabel('log2(FPKM)')
