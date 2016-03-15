%% Assessing within group correlation
% Run PCA
[coeff, score, latent] = pca(melFPKM);

figure;
stairs(cumsum(latent)/sum(latent));
ylabel('% Variance Described');
xlabel('PC');
title('Variance Described By Principal Components');

% 2 PCs can account for 95% of the variability.

% Plot first two PCs

groundtruth=[1;1;1;1;1;1;2;2;2;2;2;2]; %generating groups
group1=find(groundtruth==1);
group2=find(groundtruth==2);

figure;
plot(score(group1, 1), score(group1, 2), 'ro');
hold on
plot(score(group2, 1), score(group2, 2), 'bo');
xlabel('pc1');
ylabel('pc2');
title('Melanophore Data PCA With 2 PCs');
legend('Ablated','Unablated')
axis equal;

%% LDA
G=groundtruth';
myPCs = score(:,1:11);
N=100;

test_frac = 0.2;

for i=1:N;

 permuted = randperm(numel(myPCs(:,1)));

 test_set = permuted(1:floor(numel(myPCs(:,1))*test_frac));

 train_set = permuted(ceil((numel(myPCs(:,1))*test_frac)):end);

 LDA = fitcdiscr(myPCs(train_set,:), G(:,train_set));

 pX = predict(LDA, myPCs(test_set,:));

 lCVA = sum(G(:,test_set)-pX' ==0)/length(G(:,test_set));

 cvas(i)=lCVA;

end;

LDA = fitcdiscr(score(:,1:3), groundtruth');
pX = predict(LDA, score(:,1:3));

lCVA = sum(groundtruth'-pX' ==0)/length(groundtruth);


figure;
scatter3(score(group1, 1), score(group1, 2), score(group1, 3), 'ro');
legend('Ablated');
hold on
scatter3(score(group2, 1), score(group2, 2), score(group2, 3), 'bo');

xlabel('pc1');
ylabel('pc2');
zlabel('pc3');
title('Melanophore Data with 3 PCs');
legend('Unablated')
axis equal;

%%  QDA 

N=100;
cvas2=zeros(100);
test_frac = 0.2;


for i=1:N;

 permuted = randperm(numel(myPCs(:,1)));

 test_set = permuted(1:floor(numel(myPCs(:,1))*test_frac));

 train_set = permuted(ceil((numel(myPCs(:,1))*test_frac)):end);

 QDA = fitcdiscr(myPCs(train_set,:), G(:,train_set),...
     'DiscrimType', 'quadratic', 'Gamma', 1);

 pX2 = predict(QDA, myPCs(test_set,:));

 qCVA = sum(G(:,test_set)-pX2' ==0)/length(G(:,test_set));

 cvas2(i)=qCVA;

end;

myCVA = mean(cvas2);


%% kmeans
melkmeans=kmeans(melFPKM, 2);
melkmeans3=kmeans(melFPKM, 3);
