%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
addpath(genpath('utils'));
%I = im2double(imread('TestImage.jpg'));
%I = im2double(imread('cat.jpg'));
%I = im2double(imread('flower.jpg'));
I = im2double(imread('aman.png'));
if numel(size(I)) == 3, I = rgb2gray(I); end;
neighbors = initNeighbors(I);
%%%%%%%%%%
% parameters
%%%%%%%%%%
K = .897; 
uParams.sig = 30/255;
% for test image
%uParams.mu_b = 128; uParams.mu_f1 = 30; uParams.mu_f2 = 225;    
% for cat
%uParams.mu_b = 30; uParams.mu_f1 = 150; uParams.mu_f2 = 225;    
% for flower
%uParams.mu_b = 30; uParams.mu_f1 = 120; uParams.mu_f2 = 225;    
% for man
uParams.mu_b = 200; uParams.mu_f1 = 30; uParams.mu_f2 = 120;    
% scale [0 1] for graphcuts
uParams.mu_b = uParams.mu_b/255; uParams.mu_f1 = uParams.mu_f1/255;
uParams.mu_f2 =uParams.mu_f2/255;    

%% 1. Data term only
%% Compute costs
[unaryCosts, U, initLabels] = getUnaryCost(I(:), uParams); % N by 2
CRF.isCRF = 0;
CRF.K = K;
pairwiseCosts = getPairwiseCost(I, neighbors, CRF); % N by N
%% 2. with Pairwise term and ICM
[Uicm, icmLabels] = icm(unaryCosts, pairwiseCosts, initLabels, neighbors);

%% 3. with Simulated Annealing
[Usa, saLabels] = sim_annealing(unaryCosts, pairwiseCosts, ...
                                   initLabels, neighbors);
%% 4. CRF
sigma = 1.5;
I2 = gaussfilt(I, sigma);
[Ix Iy] = derivative5(I2, 'x', 'y');
CRF.isCRF = 1;
CRF.Mx = max(max(abs(Ix)));
CRF.My = max(max(abs(Iy)));
CRF.Ix = Ix; CRF.Iy = Iy;
CRF.Ix2 = Ix./CRF.Mx; CRF.Iy2 = Iy./CRF.My;

crfCosts = getPairwiseCost(I, neighbors, CRF);
[Ucrf, crfLabels] = icm(unaryCosts, crfCosts, saLabels, neighbors);

%% 4. CRF with graph cut algorithm 
[Ugc_a, gcLabels_a, Ugc_ab, gcLabels_ab] = graphcut(unaryCosts, crfCosts);

%%% Plot
U0 = getAllEnergy(unaryCosts, pairwiseCosts, initLabels, neighbors);
img1 = reshape(initLabels, size(I));
img2 = reshape(icmLabels, size(I));
img3 = reshape(saLabels, size(I));
img4 = reshape(crfLabels, size(I));
img5 = reshape(gcLabels_a, size(I));
img6 = reshape(gcLabels_ab, size(I));

sfigure; 
subplot(231); imagesc(img1); colormap(gray); axis image
title(sprintf('unary only E=%g', U0), 'FontSize', 15);
subplot(232); imagesc(img2); colormap(gray); axis image
title(sprintf('unary+pair ICM E=%g', Uicm), 'FontSize', 15);
subplot(233); imagesc(img3); colormap(gray); axis image
title(sprintf('unary+pair SA E=%g', Usa),'FontSize', 15);
subplot(234); imagesc(img4); colormap(gray); axis image
title(sprintf('CRF E=%g', Ucrf),'FontSize', 15);
subplot(235); imagesc(img5); colormap(gray); axis image
title(sprintf('graph cuts alpha-exp\n E=%g', Ugc_a),'FontSize', 15);
subplot(236); imagesc(img6); colormap(gray); axis image
title(sprintf('graph cuts alpha-beta swap\n E=%g', Ugc_ab),'FontSize', 15);


sfigure; imagesc(I); colormap(gray); axis image
title(sprintf('original image'),'FontSize', 15);

