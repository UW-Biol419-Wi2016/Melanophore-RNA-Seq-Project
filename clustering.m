%% Assessing within group correlation

% Ran PCA on the data

[coeff, score, latent] = pca(melFPKM);

figure;
stairs(cumsum(latent)/sum(latent));

% 95% of data is explained by 2 PCs

%% Plot the data using the first two PCs

% define ablated and unablated ground truth
groundtruth=[1;1;1;1;1;1;2;2;2;2;2;2];
group1=find(groundtruth==1);
group2=find(groundtruth==2);

% Plot
figure;
plot(score(group1, 1), score(group1, 2), 'ro');
hold on
plot(score(group2, 1), score(group2, 2), 'bo');
xlabel('pc1');
ylabel('pc2');
title('Melanophore Data PCA With 2 Groups');
axis equal;

%%