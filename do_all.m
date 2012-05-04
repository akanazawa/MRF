%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
addpath(genpath('utils'));
I = double(imread('TestImage.jpg'));
%I = double(imread('flower.jpg'));
sites = initSites(I);
%%%%%%%%%%
% parameters
%%%%%%%%%%
K = .897; 
beta = 1;

%% 1. Data term only
%% Compute costs
[unaryCosts, U, initLabels] = getUnaryCost(I(:)); % N by 2
pairwiseCosts = getPairwiseCost(I, size(I), K, []); % N by N
%% 2. with Pairwise term and ICM
[Uicm, icmLabels] = icm2(unaryCosts, pairwiseCosts, initLabels);

%% 3. with Simulated Annealing
%[Usa, saLabels] = sim_annealing(sites, initLabels, K, beta, []);
[Usa2, saLabels2] = sim_annealing2(unaryCosts, pairwiseCosts, initLabels);
%% 4. CRF
sigma = 1.5;
I2 = gaussfilt(I, sigma);
[Ix Iy] = derivative5(I2, 'x', 'y');        
CRF.K = K;
CRF.Ix = Ix; CRF.Iy = Iy;
CRF.Mx = max(max(abs(Ix)));
CRF.My = max(max(abs(Iy)));

crfCosts = getPairwiseCost(I, size(I), K, CRF);
[Ucrf2, crfLabels2] = icm2(unaryCosts, crfCosts, initLabels);

%% 4. CRF with graph cut algorithm 
[Ugc, gcLabels] = graphcut(sites, I(:), CRF);

% beta term in http://yuwing.kaist.ac.kr/courses/cs770/reading/grabcut.pdf
% CRF.meanDiffSq = 1./(mean(sum(bsxfun(@minus, I(:), I(:)').^2, 2)));
% gcLabels2 = graphcut(sites, I(:), CRF);

%%% Plot
img1 = reshape(initLabels, size(I));
img2 = reshape(icmLabels, size(I));
img3 = reshape(saLabels, size(I));
img4 = reshape(crfLabels, size(I));
img5 = reshape(gcLabels, size(I));
% img6 = reshape(gcLabels2, size(I));

sfigure; 
subplot(231); imagesc(img1); colormap(gray); axis image
title('unary only');
subplot(232); imagesc(img2); colormap(gray); axis image
title('unary+pair ICM');
subplot(233); imagesc(img3); colormap(gray); axis image
title('unary+pair SA');
subplot(234); imagesc(img4); colormap(gray); axis image
title('CRF');
subplot(235); imagesc(img5); colormap(gray); axis image
title('CRF-graph cuts');
subplot(236); imagesc(img6); colormap(gray); axis image
title('CRF-graph cuts with different pairwise potential');

