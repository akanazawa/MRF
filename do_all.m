%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
addpath(genpath('utils'));
%I = double(imread('TestImage.jpg'));
I = double(imread('sil.jpg'));
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
[Uicm, icmLabels] = icm2(unaryCosts, pairwiseCosts, initLabels, sites);

%% 3. with Simulated Annealing
[Usa, saLabels] = sim_annealing2(unaryCosts, pairwiseCosts, ...
                                   initLabels, sites);
%% 4. CRF
sigma = 1.5;
I2 = gaussfilt(I, sigma);
[Ix Iy] = derivative5(I2, 'x', 'y');        
CRF.K = K;
CRF.Ix = Ix; CRF.Iy = Iy;
CRF.Mx = max(max(abs(Ix)));
CRF.My = max(max(abs(Iy)));
crfCosts = getPairwiseCost(I, size(I), K, CRF);
[Ucrf, crfLabels] = icm2(unaryCosts, crfCosts, initLabels, sites);

%% 4. CRF with graph cut algorithm 
[Ugc, gcLabels] = graphcut(unaryCosts, crfCosts);

% beta term in http://yuwing.kaist.ac.kr/courses/cs770/reading/grabcut.pdf
% CRF.meanDiffSq = 1./(mean(sum(bsxfun(@minus, I(:), I(:)').^2, 2)));
% gcLabels2 = graphcut(sites, I(:), CRF);

%%% Plot
U0 = getAllEnergy(unaryCosts, pairwiseCosts, initLabels);
img1 = reshape(initLabels, size(I));
img2 = reshape(icmLabels, size(I));
img3 = reshape(saLabels, size(I));
img4 = reshape(crfLabels, size(I));
img5 = reshape(gcLabels, size(I));
% img6 = reshape(gcLabels2, size(I));

sfigure; 
subplot(231); imagesc(img1); colormap(gray); axis image
title(sprintf('unary only E=%g', U0));
subplot(232); imagesc(img2); colormap(gray); axis image
title(sprintf('unary+pair ICM E=%g', Uicm));
subplot(233); imagesc(img3); colormap(gray); axis image
title(sprintf('unary+pair SA E=%g', Usa));
subplot(234); imagesc(img4); colormap(gray); axis image
title(sprintf('CRF E=%g', Ucrf));
subplot(235); imagesc(img5); colormap(gray); axis image
title(sprintf('CRF-graph cuts E=%g', Ugc));
subplot(236); imagesc(I); colormap(gray); axis image
title(sprintf('original image'));

