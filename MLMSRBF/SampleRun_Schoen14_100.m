% sample run of the MLMSRBF algorithm

% set random seed
rng(23,'twister');

% problem information
probinfo.funcname = 'Schoen14_100';
probinfo.lb = zeros(1,14);
probinfo.ub = ones(1,14);

% specify initial point or generate uniformly at random throughout search space
probinfo.userinitpoints = probinfo.lb + rand(1,14).*(probinfo.ub - probinfo.lb);

% algorithm settings
algparams.rbftype = 'thinplate';  % setting used in the paper but may be replaced by 'cubic'
algparams.maxeval = 500;

% run algorithm
outputinfo = RunMLMSRBF(probinfo, algparams);