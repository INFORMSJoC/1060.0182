% MLMSRBF.m implements the Multistart Local Metric Stochastic RBF
% (MLMSRBF) method for finding the global minimum of a function defined
% over a box. Described in the paper:
%
% R.G. Regis, C.A. Shoemaker. A stochastic radial basis function method 
% for the global optimization of expensive functions. INFORMS Journal on
% Computing, Vol. 19, No. 4, pp. 497-509, 2007.
%
% Written by: Rommel G. Regis
% Last updated: 07/03/07, 06/24/10, 06/29/10
%
% Syntax:
% outputinfo = MLMSRBF(probinfo, algparams)
%
%
% Input:
%
%   probinfo.funcname: Name of function represented as a string
%   probinfo.lb: row vector of lower bounds
%   probinfo.ub: row vector of upper bounds
%   probinfo.optvalue: global minimum value of the function over the domain
%   probinfo.userinitpoints: initial points provided by the user
%
%   algparams.rbftype: type of RBF (options: 'cubic', 'thinplate', 
%                                   'multiquadric', 'Gaussian')
%
%   algparams.maxeval: maximum number of function evaluations 
%   (The algorithm terminates when the total number of calls to the
%   objective function exceeds maxeval.)
%
%   algparams.numdesignpts: number of experimental design points for every
%   run of LMSRBF (must be >= dim+1)
%
%   algparams.accuracy: required percent accuracy for the best function
%   value (The algorithm terminates when the best function value found
%   has a relative error that is less than accuracy. Recommended Value = 1)
%
%   algparams.transform: indicates whether the high function values will be
%   replaced by the median of all available function values before fitting
%   the RBF model (Set to 'T' (replace high function values) or 'N' (use
%   original function values. Recommended Value = 'T'.)
%
%
% Output:
%
%   outputinfo.bestsolution: best solution vector found
%
%   outputinfo.bestvalue: function value of best solution vector (i.e. the
%   best function value found
%
%   outputinfo.numevals: actual number of calls to the objective function
%
%   outputinfo.proctime: total elapsed time (in seconds)
%
%   outputinfo.totalfevaltime: total time (in seconds) spent on all function
%   evaluations
% 
%   outputinfo.success: boolean-valued variable indicating whether the best
%   function value found has the required accuracy (1 if true, 0 otherwise).
%
%   outputinfo.history: history matrix that summarizes the entire history 
%   and trajectory of the algorithm from the beginning until it terminated.
%   It is a matrix of size numevals x (d+4), where d is the dimension of
%   the problem. Each row of the history matrix corresponds to one function
%   evaluation. In particular, for the nth function evaluation, we have:
%   history(n,1) is the LMSRBF run number
%   history(n,2) is the function value of the nth function evaluation point
%   history(n,3:d+2) is the nth function evaluation point
%   history(n,d+3) is the time (in seconds) to execute the nth function
%   evaluation
%   history(n,d+4) is time elapsed (in seconds) up to the nth function
%   evaluation excluding the time spent on all function evaluations so far
%   (i.e., excluding the total time for the first n function evaluations).

function outputinfo = MLMSRBF(probinfo, algparams)

% gather the input
funcname = probinfo.funcname;
xmin = probinfo.lb;
xmax = probinfo.ub;
Domain = [xmin; xmax];
dim = size(Domain,2);       % number of continuous variables
if (isfield(probinfo,'optvalue') == 0)
    optvalue = -1e300;      % lower bound for optimal value
else
    optvalue = probinfo.optvalue;
end
if (isfield(probinfo,'userinitpoints') == 0)
    UserInitPoints = [];
else
    disp('userinitpts provided');
    UserInitPoints = probinfo.userinitpoints;
end

if (isfield(algparams,'rbftype') == 0)
    phifunction = 'thinplate';      % setting used in original paper
else
    phifunction = algparams.rbftype;
end
parameter = [];
polynomial = DeterminePolynomialTail(phifunction);
maxeval = algparams.maxeval;
if (isfield(algparams,'numdesignpts') == 0)
    numdesignpts = 2*(dim + 1);     % setting used in original paper
else
    numdesignpts = max(algparams.numdesignpts, dim+1);  % ensure there is enough design points
end
if (isfield(algparams,'accuracy') == 0)
    accuracy = 0;       % setting used in original paper (allows algorithm to exhaust maxeval)
else
    accuracy = algparams.accuracy;
end
if (isfield(algparams,'transform') == 0)
   transform = 'T';     % setting used in original paper
else
   transform = algparams.transform;
end

% set the timer
time0 = clock;

% target value for the given accuracy
targetvalue = optvalue + abs(optvalue)*(accuracy/100);

% number of starts
numstart = 0;

bestsolution = [];
bestvalue = inf;
numevals = 0;
totalfevaltime = 0;
history = [];

while ((bestvalue > targetvalue) && (numevals < maxeval))
    
    numstart = numstart + 1;
    fprintf('\n\nLMSRBF run number: %d\n\n',numstart);
    
    % select space-filling experimental design
    if (numstart == 1)
        DesignPoints = GetInitEvalPoints(UserInitPoints, numdesignpts, xmin, xmax, polynomial, 'none');
    else
        disp('Generating a new symmetric Latin hypercube design (SLHD).');
        TempDesignPoints = SelectUnisolventSLHD(dim, numdesignpts, polynomial);
        DesignPoints = Rescale(TempDesignPoints, xmin, xmax);
    end
    
    % run LMSRBF
    [bestsolution_temp,bestvalue_temp,numevals_temp,proctime_temp,totalfevaltime_temp,success_temp,history_temp] = ...
        LMSRBF(funcname,Domain,phifunction,parameter,polynomial,DesignPoints,optvalue,accuracy,maxeval-numevals,transform);
        
    % update the history matrix
    history_temp(:,1) = numstart*ones(numevals_temp,1);
    if (numstart == 1)
        history = history_temp;
    else
        history = [history; history_temp];
        history((numevals+1):(numevals+numevals_temp),dim+4) = history(numevals,dim+4) + history_temp(:,dim+4);
    end
        
    % update the results
    if (bestvalue_temp < bestvalue)
        bestsolution = bestsolution_temp;
        bestvalue = bestvalue_temp;
    end    
    numevals = numevals + numevals_temp;
    totalfevaltime = totalfevaltime + totalfevaltime_temp;
    success = success_temp;
     
end

fprintf('\n\nTotal number of LMSRBF runs: %d\n\n',numstart);
proctime = etime(clock,time0);

% gather the output
outputinfo.bestsolution = bestsolution;
outputinfo.bestvalue = bestvalue;
outputinfo.numevals = numevals;
% outputinfo.proctime = proctime;
% outputinfo.totalfevaltime = totalfevaltime;
outputinfo.history = history;
% outputinfo.success = success;





% subfunctions

% determines the corresponding polynomial tail for a particular radial function

function polynomial = DeterminePolynomialTail(phifunction)

switch (phifunction)
    case {'cubic','thinplate'}
        polynomial = 'linear';
    case 'multiquadric'
        polynomial = 'constant';
    case 'Gaussian'
        polynomial = 'none';
    otherwise
        disp('ERROR: Unknown RBF type!');
        return;
end



% Notes: numexpdesign >= dim + 1
%        UserInitPoints should already be in the right scale
%        designcriterion: 'none' (setting used in original paper) or 'maximin'

function InitPoints = GetInitEvalPoints(UserInitPoints, numdesignpts, xmin, xmax, polynomial, designcriterion)

if (IsUnisolventLinear(UserInitPoints) == 1)
    % if points provided by the user are suitable for RBF with linear tail
    % interpolation then use them
    InitPoints = UserInitPoints;
    disp('Points provided by the user are suitable initial points and will be used by the algorithm.');
else
    % create a unisolvent maximin SLHD
    dim = length(xmin);
    if (strcmp(designcriterion,'maximin') == 1)
        disp('Generating an approximately optimal symmetric Latin hypercube design (SLHD).');
        maxrandSLHD = min(1000, 10*numdesignpts*dim);
        TempDesignPoints = SelectUnisolventMaximinSLHD(dim, numdesignpts, polynomial, maxrandSLHD);
    else
        disp('Generating a symmetric Latin hypercube design (SLHD)');
        TempDesignPoints = SelectUnisolventSLHD(dim, numdesignpts, polynomial);
    end
    DesignPoints = Rescale(TempDesignPoints, xmin, xmax);
        
    % combine user points with design points
    CombPoints = [UserInitPoints; DesignPoints];
    
    % remove any repeated rows
    minxrange = min(xmax - xmin);
    tolerance = 0.0005*minxrange*sqrt(dim);
    % [InitPoints, numdistinct] = RemoveRepeatedRows(CombPoints, tolerance);
    [InitPoints, ~] = RemoveRepeatedRows(CombPoints, tolerance);
    
end



% checks if the given set of points is suitable for RBF-linear interpolation
% Note: The points are given as row vectors.

function booloutput = IsUnisolventLinear(Points)

[numpts, dim] = size(Points);
P = [ones(numpts,1), Points];
if (rank(P) == dim+1)
    booloutput = 1;
else
    booloutput = 0;
end



% selects a unisolvent maximin SLHD

function DesignPoints = SelectUnisolventMaximinSLHD(dim, numdesignpts, polynomial, maxrandSLHD)

if (strcmp(polynomial,'linear') == 1)
    P = [0];
    while (rank(P) ~= dim+1)
        DesignPoints = optimizeSLHD(dim, numdesignpts, maxrandSLHD);
        P = [ones(numdesignpts,1), DesignPoints];
    end
else
    DesignPoints = optimizeSLHD(dim, numdesignpts, maxrandSLHD);
end



% optimizeSLHD creates an approximately optimal symmetric latin hypercube design.
% d is the dimension of the input and m is the number of initial points to be selected.

function [InitialPoints,Permutations,BestPairwiseDist] = optimizeSLHD(d,m,numtrials)

numdist = m*(m-1)/2;
BestPairwiseDist = zeros(1,numdist);
SortedBestPairwiseDist = zeros(1,numdist);

numreplace = 0;
for k = 1:numtrials
    [I_current,P_current] = SLHDstandard(d,m);
    CurrentPairwiseDist = ComputePairwiseDistanceVector(I_current);
    SortedCurrentPairwiseDist = sort(CurrentPairwiseDist);
    if (isGreater(SortedCurrentPairwiseDist,SortedBestPairwiseDist) == 1)
        BestPairwiseDist = CurrentPairwiseDist;
        SortedBestPairwiseDist = SortedCurrentPairwiseDist;
        I_best = I_current;
        P_best = P_current;
        numreplace = numreplace + 1;
    end
end
InitialPoints = I_best;
Permutations = P_best;
fprintf('\nNumber of times experimental design was replaced: %d\n',numreplace);



% calculates the pairwise distances between points
% same output format as the function pdist from Matlab's Statistics Toolbox

function pdistvec = ComputePairwiseDistanceVector(DataPoints)

numpoints = size(DataPoints, 1);
numpairs = numpoints*(numpoints - 1)/2;
pdistvec = zeros(1,numpairs);
pos = 1;
for i = 1:(numpoints-1)
    xtemp = ones(numpoints-i, 1)*DataPoints(i, :);
    tempvec = sqrt(sum((xtemp-DataPoints(i+1:numpoints, :)).^2, 2));
    pdistvec(pos:(pos+numpoints-i-1)) = tempvec';
    pos = pos+numpoints-i;
end



% compares two vectors x and y

function output = isGreater(x,y)

if (length(x) ~= length(y))
    disp('Error: Input vectors must have the same length.');
    return;
end

n = length(x);
if (n == 0)
    disp('Error: Input vectors are empty.');
    return;
end

output = 2;     % 0 = no, 1 = yes, 2 = undecided
pos = 1;

while ((output == 2) && (pos <= n))
    if (x(pos) > y(pos))
        output = 1;
    else
        if (x(pos) < y(pos))
            output = 0;
        else
            pos = pos + 1;
        end
    end
end

if (output == 2)    % results when the two vectors are equal
    output = 0;
end



% selects a unisolvent SLHD

function DesignPoints = SelectUnisolventSLHD(d, m, polynomial)

if (strcmp(polynomial,'linear') == 1)
    rank_P = 0;
    while (rank_P ~= d+1)
        DesignPoints = SLHDstandard(d, m);
        P = [ones(m,1), DesignPoints];
        rank_P = rank(P);
    end
else
    DesignPoints = SLHDstandard(d, m);
end



% SLHD creates a symmetric latin hypercube design. d is the dimension of the input and
% m is the number of initial points to be selected.

function [InitialPoints,P] = SLHDstandard(d,m)

delta = (1/m)*ones(1,d);

X = zeros(m,d);
for j = 1:d
    for i = 1:m
        X(i,j) = ((2*i-1)/2)*delta(j);
    end
end

P = zeros(m,d);
P(:,1) = (1:m)';
if (mod(m,2) == 0)
   k = m/2;
else
   k = (m-1)/2;
   P(k+1,:) = (k+1)*ones(1,d);
end

for j = 2:d
   P(1:k,j) = randperm(k)';
   for i = 1:k
      if (rand(1) <= 0.5)
         P(m+1-i,j) = m+1-P(i,j);
      else
         P(m+1-i,j) = P(i,j);
         P(i,j) = m+1-P(i,j);
      end
   end
end

InitialPoints = zeros(m,d);
for j = 1:d
    for i = 1:m
        InitialPoints(i,j) = X(P(i,j),j);
    end
end



% adjust points to the correct scale

function RescaledPoints = Rescale(Points, xmin, xmax)

[numpoints, dim] = size(Points);
RescaledPoints = zeros(numpoints, dim);
xrange = xmax - xmin;
for i = 1:numpoints
    RescaledPoints(i,:) = xmin + xrange.*Points(i,:);
end



% removes repeated rows in InputMatrix

function [OutputMatrix, numdistinct] = RemoveRepeatedRows(InputMatrix, tolerance)

[numrows, numcols] = size(InputMatrix);
OutputMatrix = zeros(numrows, numcols);
numdistinct = 0;
row = 1;
while (row <= numrows)
    while ( isDiff(InputMatrix(row,:),OutputMatrix(1:numdistinct,:),tolerance) == 0 )
        row = row + 1;
        if (row > numrows)
            break;
        end
    end
    if (row > numrows)
        break;
    end
    numdistinct = numdistinct + 1;
    OutputMatrix(numdistinct,:) = InputMatrix(row,:);
    row = row + 1;
end
OutputMatrix = OutputMatrix(1:numdistinct,:);



% checks if the row vector x is different from each row of the matrix A

function output = isDiff(x,A,tolerance)

output = 1;
for i = 1:size(A,1)
    if (norm(x-A(i,:)) < tolerance)
        output = 0;
        break;
    end
end