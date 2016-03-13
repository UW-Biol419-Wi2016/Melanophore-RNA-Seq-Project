%% Start of the mathworks workflow
load melRNAseq.mat;
melData(1:10,:)

% assign groups for ablated and unablated (our 2 conditions)
ablated=melFPKM(1:6,:);
unablated=melFPKM(7:12,:);

% Workflow suggests that we can make a summary table of reads that were
% assigned vs unassigned if our files came with that information.

% Chromosome mapping portion is dependent on having chromosome position
% information in the starting data set... Shoot.

% Normalization
% estimate pseudo-reference with geometric mean row by row
melFPKMalt1=melFPKM';
pseudoRefSample = geomean(melFPKMalt1,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,melFPKMalt1(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);

%figure;

%subplot(2,1,1)
%maboxplot(log2(melFPKMalt),'title','Raw Read Count','orientation','horizontal')
%ylabel('sample')
%xlabel('log2(FPKM)')

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
mF4=melFPKM';
mF2 = melFPKMalt'; % it didn't work with the non-transposed version

% find mean (with transpose)
meanUAnew = mean(mF2(:, 7:12),2);
meanAnew = mean(mF2(:, 1:6),2);

meanUA = mean(mF4(:, 7:12),2);
meanA = mean(mF4(:, 1:6),2);

% convert zeros to NaN
% meanAnew(meanA ==0) = NaN;
% meanUAnew(meanUA ==0) = NaN;

% find dispersion (with transpose)
dispUA = std(mF2(:,7:12),0,2) ./ meanUAnew;
dispA = std(mF2(:,1:6),0,2) ./ meanAnew;

% plot on a log-log scale (with transpose)
figure;
loglog(meanUAnew, dispUA, 'or');
hold on;
loglog(meanAnew, dispA, 'ob');
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

%p=polyfit(log2(meanUA),log2(meanA),1);
%y=polyval(p,log2(meanA));
x=(1:20);
y=x-1;
>>>>>>> master
% scatter plot of means
figure;
hold on;
plot(x,y, 'r-');
plot(log2(meanUA),log2(meanA), 'o');

%p=polyfit(log2(meanUA),log2(meanA),1);
%y=polyval(p,log2(meanA));
x=(1:20);
y=x-1;
% scatter plot of means
figure;
hold on;
plot(x,y, 'r-');
plot(log2(meanUAnew),log2(meanAnew), 'o');
plot(log2(meanUAnew),log2(meanAnew), 'o');
xlabel('UnAblated Ctrl');
ylabel('Ablated');

%% Fold Change

% compute the mean and the log2FoldChange
meanBase = (meanUAnew + meanAnew) / 2;
foldChange = meanAnew ./ meanUAnew;
log2FC = log2(foldChange);

% plot mean vs. fold change (MA plot)
mairplot(meanAnew, meanUAnew,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')

