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

%% To plot the mean FPKM values between the two groups

mF2 = melFPKM'; % it didn't work with the non-transposed version

% find mean (with transpose)
meanUA = mean(mF2(:, 7:12),2);
meanA = mean(mF2(:, 1:6),2);

% find dispersion (with transpose)
dispUA = std(mF2(:,7:12),0,2) ./ meanUA;
dispA = std(mF2(:,1:6),0,2) ./ meanA;

% plot on a log-log scale (with transpose)
figure;
loglog(meanUA, dispUA, '.r');
hold on;
loglog(meanA, dispA, '.b');
xlabel('log2(mean)');
ylabel('log2(Dispersion)');
legend('UnAblated', 'Ablated', 'Location', 'southwest');
