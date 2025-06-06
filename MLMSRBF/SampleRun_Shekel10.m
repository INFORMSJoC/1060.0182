% sample run of the MLMSRBF algorithm

% set random seed
rng(23,'twister');

% problem information
probinfo.funcname = 'Shekel10';
probinfo.lb = zeros(1,4);
probinfo.ub = 10*ones(1,4);

% specify initial point or generate uniformly at random throughout search space
probinfo.userinitpoints = [5.2, 3.0, 1.9, 2.6; 
                           1.2, 5.7, 3.4, 8.6];
% probinfo.userinitpoints = probinfo.lb + rand(1,4).*(probinfo.ub - probinfo.lb);

% algorithm settings
algparams.rbftype = 'thinplate';  % setting used in the paper but may be replaced by 'cubic'
algparams.maxeval = 200;

% run algorithm
outputinfo = RunMLMSRBF(probinfo, algparams);