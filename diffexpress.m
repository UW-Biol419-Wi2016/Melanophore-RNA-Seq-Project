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
%% Removing FPKM values below 1

melFPKMalt=melFPKM;
for i=1:38125;
    for j=1:12;
        if melFPKMalt(j,i)<1;
            melFPKMalt(j,i)=NaN;
        end;
    end;
end;
%% To plot the mean FPKM values between the two groups

mF2 = melFPKMalt'; % it didn't work with the non-transposed version

% find mean (with transpose)
meanUA = mean(mF2(:, 7:12),2);
meanA = mean(mF2(:, 1:6),2);

% find dispersion (with transpose)
dispUA = std(mF2(:,7:12),0,2) ./ meanUA;
dispA = std(mF2(:,1:6),0,2) ./ meanA;

% plot on a log-log scale (with transpose)
figure;
loglog(meanUA, dispUA, 'or');
hold on;
loglog(meanA, dispA, 'ob');
xlabel('log2(mean)');
ylabel('log2(Dispersion)');
legend('UnAblated', 'Ablated', 'Location', 'southwest');


% convert zeros to NaN
meanA(meanA ==0) = NaN;
meanUA(meanUA ==0) = NaN;

mdl=fitlm(meanUA,meanA);

% linear regression
[p,ErrorEst] = polyfit(log2(meanUA),log2(meanA),1);
pop_fit = polyval(p, log2(meanUA), ErrorEst);


% scatter plot of means
figure;
plot(log2(meanUA), pop_fit, '-');
hold on;
plot(log2(meanUA), log2(meanA),'o');
% plot(log2(meanUA),log2(meanA), 'o');
xlabel('UnAblated Ctrl');
ylabel('Ablated');

x=(1:20);
y=x-0.5;

% scatter plot of means
figure;
hold on;
plot(x,y, 'r-');
plot(log2(meanUA),log2(meanA), 'o');
xlabel('UnAblated Ctrl');
ylabel('Ablated');

% 

%% Fold Change

% compute the mean and the log2FoldChange
meanBase = (meanUA + meanA) / 2;
foldChange = meanA ./ meanUA;
log2FC = log2(foldChange);

% plot mean vs. fold change (MA plot)
mairplot(meanA, meanUA,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')

%% create table with gene statistics 

geneTable = table(meanBase, meanA, meanUA, foldChange, log2FC);
geneTable.properties.RowNames = zgenes.tracking_id