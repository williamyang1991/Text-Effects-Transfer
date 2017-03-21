function [optS] = init_opt()

% INIT_OPT
% Initialize parameters for analysis and synthesis

fprintf('- Initialize parameters \n');

optS = [];

% === enable randomness ===

optS.lockRandOn = false;            
optS.lockRandSeed = 0;

% === optimal patch scale detection ===

optS.histNum  = 100;                % quantify all distances into 100 bins
optS.numPSizeLvl = 5;               % the max scale L = 5

% === Method configuration ===
optS.useBiasCorrection = true;      % permit patch color correction
optS.initByDistOn = true;           % initialize NNF with Distribution Term
optS.visOn = false;                 % show intermediate result during the process

% === Weighting parameters ===

optS.lambdaDist  = 0.1;             % ¦Ë1 = 0.1
optS.lambdaRep   = 0.005;           % ¦Ë2 = 0.005
optS.lambdaText  = 10;              % ¦Ë3 = 10

% === Patch size ===
optS.pSize = 5;                        % Patch size (odd number), use larger patch for more coherent region
optS.pRad  = floor(optS.pSize/2);      % Patch radius
optS.pNumPix = optS.pSize*optS.pSize;  % Number of pixels in a patch
optS.pMidPix = round(optS.pNumPix/2);  % The center of the patch

% === Multi-resolution parameters ===
optS.numPyrLvl = 10;                   % Number of coarse to fine layer
optS.coarestImgSize = 32;              % The size of the smallest image in the pyramid
optS.resampleKernel = 'lanczos3';      % Resampling kernel: 'lanczos3', 'bicubic', 'bilinear'
optS.useLogScale = 1;                  % Use log scales or linear scales for downsampling

% Weighting parameters for patch match, larger weight put more emphasis on
% pixels near to text
optS.wDist = 2.^linspace(1, 1, optS.numPyrLvl);

optS.topLevel = 1;                     % Which level to stop
optS.propOnlyLevel = 0;

% === Parameters for patch transformation ===

optS.maxBias =  0.05;                   % Maximum bias compensation
optS.minBias = -optS.maxBias;           % Mininum bias compensation
optS.angleRad = 0.5;                    % Maximum patch rotation

% === Number of iterations per level ===
optS.numIter    = 30;                  % The initial iteration number, large hole might require
% more iterations
optS.numIterDec = optS.numIter/optS.numPyrLvl;   % Number of decrements
optS.numIterMin = optS.numIter/optS.numPyrLvl;   % Minimum number of iterations
optS.numPassPerIter = 1;

% === Precomputed patch position in the reference position ===
[X, Y] = meshgrid(-optS.pRad:optS.pRad, -optS.pRad:optS.pRad);
optS.refPatchPos = single(cat(1, X(:)', Y(:)', ones(1, optS.pSize*optS.pSize)));

% === Propagation directions ===
optS.propDir = [1 0; 0 1; -1 0; 0 -1]';

optS.rsThres = 0.05;                   % Random sampling threshold for early termination
optS.voteUpdateW = 20;                 % Smaller number for quick voting update

optS.numRegSample  = 1;                % Number of regularity-guided sampling per iteration
optS.numRandSample = 5;                % Number of coarse-to-fine random sampling per iteration

% [To-Do] robust norm, e.g., huber
optS.costType = 'L2';                  % Patch appearance cost, other option: 'L2'

end