% Local Metric Stochastic Radial Basis Function Method
% Note: InitialPoints should already be in the right scale

function [solution,value,numevals,proctime,totalfevaltime,success,history] ...
    = LMSRBF(funcname,Domain,phifunction,parameter,polynomial,...
    InitialPoints,optvalue,accuracy,maxeval,transform)

% set problem parameters
global d xmin xmax xrange minxrange targetvalue;                   
d = size(Domain,2);                                     % number of continuous variables
xmin = Domain(1,:);                                     % lower bounds on the variables
xmax = Domain(2,:);                                     % upper bounds on the variables
xrange = xmax - xmin;                                   % ranges of the variables
minxrange = min(xrange);                                % side length of the hypercube
targetvalue = optvalue + abs(optvalue)*(accuracy/100);  % target value for the given accuracy

% set RBF parameters
global rbftype gamma polytail;             
rbftype = phifunction;              % type of rbf
gamma = parameter*minxrange;        % parameter for the rbf part
polytail = polynomial;              % type of polynomial tail

% set tolerance parameters
global tolerance maxconstrviol;
tolerance = 0.0005*minxrange*sqrt(d);
maxconstrviol = (1e-6)*minxrange;

% initialize variables, vectors and matrices
global D F fevaltime elapsedtime n;
D = zeros(maxeval, d);              % stores the evaluation points
F = zeros(maxeval, 1);              % stores the function values of the evaluated points
fevaltime = zeros(maxeval, 1);      % time for function evaluations
elapsedtime = zeros(maxeval, 1);    % time elapsed so far
n = 0;                              % counter for the number of function evaluations

% initialize the current best solution
global xbest Fbest;
xbest = [];
Fbest = inf;

% determine the number of initial evaluation points
global numinitpts;
numinitpts = size(InitialPoints,1);


%%%%%%% start of the algorithm %%%%%%

% set the timer
global time0;
time0 = clock;

% evaluate the experimental design points
numinitpts = min(numinitpts, maxeval);
for i = 1:numinitpts
    EvaluateFunction(funcname, InitialPoints(i,:));
end

% check if maxeval has been reached
if (n >= maxeval)
    history = CreateHistoryMatrix;
    [solution, value, numevals, proctime, totalfevaltime, success] = GetOutputValues;
    return;    
end

% do leave-one-out cross-validation (LOOCV) to determine gamma for
% multiquadric or Gaussian
PerformLOOCV(numinitpts);

% generate initial RBF matrices
global PHI phi0 P pdeg pdim;
[PHI, phi0, P, pdeg, pdim] = ComputeInitialRBFMatrices(maxeval, numinitpts, rbftype, gamma, polytail);

% set local stochastic RBF algorithm parameters (global variables will be set here)
SetLocalStochRBFAlgorithmParams;

% initialize control parameters (global variables will be set here)
global localminflag;
InitializeLocalStochRBFControlParams;

% global variables needed by functions in the while loop
% Ftransform contains the function values of the previously evaluated points (possibly transformed)
% lambda and ctail are the coefficients of the RBF model
% Fbest_old is the previous value of Fbest
global Ftransform lambda ctail Fbest_old;

% initialize the number of iterations
iterctr = 0;                                    


% iteration steps

% while ((n < maxeval) && (localminflag == 0))
while (((n < maxeval) && (localminflag == 0)) && (Fbest > targetvalue))
   
    iterctr = iterctr + 1;
    fprintf('\nIteration %d:\n',iterctr);
        
    % replace large function values by the median of all available function values
    Ftransform = TransformFunctionValues(transform);
    
     % fit RBF approximation model
    [lambda, ctail] = FitRBF;
        
    % select function evaluation point
    [xselected, normval] = SelectLocalStochRBFPoint(iterctr);
        
    % perform function evaluation at the selected point
    Fbest_old = Fbest;    
    Fselected = EvaluateFunction(funcname, xselected);
    % fprintf('current value: %.8e\n',Fselected);
    % fprintf('best value so far: %.8e (target: %.8e)\n',Fbest,targetvalue);
    fprintf('current value: %.8e (best value so far: %.8e)\n',Fselected, Fbest);
    
    % update local stochastic RBF control parameters
    UpdateLocalStochRBFControlParams;
    
    % do the following only if planning to do another iteration
    % if ((n < maxeval) && (localminflag == 0))
    if (((n < maxeval) && (localminflag == 0)) && (Fbest > targetvalue))   
        UpdateRBFMatrices(xselected, normval);
    end
        
end

% create the history matrix
history = CreateHistoryMatrix;

% get output values
[solution, value, numevals, proctime, totalfevaltime, success] = GetOutputValues;





% subfunctions

% evaluates the function at a given point and updates the global variables
% n, D, F, fevaltime, elapsedtime, xbest and Fbest
%
% Inputs:
% (1) funcname is a string representing the name of the function
% (2) evalpoint is the given point
%
% Output:
% funcvalue is the value of the function at evalpoint

function funcvalue = EvaluateFunction(funcname, evalpoint)

global D F fevaltime elapsedtime n xbest Fbest time0;

time1 = clock;
funcvalue = feval(funcname, evalpoint);
n = n + 1;
fevaltime(n) = etime(clock, time1);
elapsedtime(n) = etime(clock, time0);
if (funcvalue < Fbest)
    xbest = evalpoint;
    Fbest = funcvalue;
end
D(n,:) = evalpoint;
F(n) = funcvalue;



% creates the history matrix for the entire run

function History = CreateHistoryMatrix

global D F fevaltime elapsedtime n;

optime = zeros(n,1);
for i = 1:n
    optime(i) = elapsedtime(i) - sum(fevaltime(1:i));
end
History = [zeros(n,1), F(1:n), D(1:n,:), fevaltime(1:n), optime];



% performs leave-one-out cross validation

function PerformLOOCV(numinitpts)

global d minxrange rbftype gamma polytail D F;

if ( ((strcmp(rbftype,'multiquadric') == 1) || (strcmp(rbftype,'Gaussian') == 1)) && (isempty(gamma) == 1) )
    if (strcmp(rbftype,'Gaussian') == 1)
        mingamma = (minxrange*sqrt(d))/(4*sqrt(log(10)));
        maxgamma = (0.005/sqrt(-log(0.9995)))*minxrange*sqrt(d);
    else
        mingamma = 0.005*sqrt(d)*minxrange;
        maxgamma = (sqrt(d)/2)*minxrange;
    end
    gammarange = mingamma:((maxgamma-mingamma)/20):maxgamma;
    [gamma,~] = LOOCV(D(1:numinitpts,:),F(1:numinitpts),rbftype,polytail,gammarange);
end



% implements LOOCV to determine gamma parameter for RBF methods

function [gammaparam, rmse] = LOOCV(D, F, rbftype, polytail, gammavalues)

numgamma = length(gammavalues);
rmse = zeros(1,numgamma);

[m, d] = size(D);

switch polytail
case 'none'
    pdeg = -1;
    pdim = 0;
    P = [];
case 'constant'
    pdeg = 0;
    pdim = 1;
    P = ones(m,1);
case 'linear'
    pdeg = 1;
    pdim = d + 1;
    P = [ones(m,1),D];
case 'quadratic'
    pdeg = 2;
    pdim = ((d+1)*(d+2))/2;
    P = [ones(m,1),D,zeros(m,d*(d+1)/2)];
    columnpos = d+1;
    for i = 1:d
        for j = i:d
            columnpos = columnpos + 1;            
            P(:,columnpos) = D(:,i).*D(:,j);
        end
    end
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end

% determine pairwise distance between points
PairwiseDistance = zeros(m);
for i = 1:(m-1)
    xtemp = ones(m-i,1)*D(i,:);
    PairwiseDistance(i,i+1:m) = sqrt(sum((xtemp-D(i+1:m,:)).^2,2))';
end


for h = 1:numgamma
    
    PHI = zeros(m);
    phi0 = phi(0,rbftype,gammavalues(h));
    for i = 1:m
        PHI(i,i) = phi0;
        PHI(i,i+1:m) = phi(PairwiseDistance(i,i+1:m),rbftype,gammavalues(h));
        PHI(i+1:m,i) = PHI(i,i+1:m)';
    end
    
    sumsquarederror = 0;
        
    for k = 1:m
            
        PHItemp = [PHI(1:k-1,:); PHI(k+1:m,:)];
        PHItemp = [PHItemp(:,1:k-1), PHItemp(:,k+1:m)];
            
        if (strcmp(polytail,'none') == 1)
            Ptemp = [];
        else
            Ptemp = [P(1:k-1,:); P(k+1:m,:)];
        end
            
        Dtemp = [D(1:k-1,:); D(k+1:m,:)];
        Ftemp = [F(1:k-1); F(k+1:m)];
            
        coefftemp = [PHItemp,Ptemp;Ptemp',zeros(pdim)]\[Ftemp;zeros(pdim,1)];
        lambdatemp = coefftemp(1:m-1);
        ctemp = coefftemp(m:m-1+pdim);
            
        sumsquarederror = sumsquarederror + (RBFvalue(D(k,:),Dtemp,rbftype,gammavalues(h),polytail,lambdatemp,ctemp) - F(k))^2;
            
    end
        
    rmse(h) = sqrt(sumsquarederror/m);
        
end
            
[~,index] = min(rmse);
gammaparam = gammavalues(index);
fprintf('gamma parameter obtained by LOOCV: %f\n\n',gammaparam);



% RBF value
%
% syntax:
% [f,g] = RBFvalue(x,D,rbftype,gamma,polytail,lambda,ctail)

function [f,g] = RBFvalue(x,D,rbftype,gamma,polytail,lambda,ctail)

% x is a row vector of length d

[n,d] = size(D);
xtemp = ones(n,1)*x;
norm_value = sqrt(sum((xtemp-D).^2,2));
phi_vector = phi(norm_value,rbftype,gamma);

switch polytail
case 'none'
    polypart = 0;
case 'constant'
    polypart = ctail;
case 'linear'
    polypart = [1,x]*ctail;
case 'quadratic'
    temp = [1,x,zeros(1,(d*(d+1))/2)];
    columnpos = d+1;
    for i = 1:d
        for j = i:d
            columnpos = columnpos + 1;
            temp(columnpos) = x(i)*x(j);
        end
    end
    polypart = temp*ctail;
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end

f = (lambda')*phi_vector + polypart;


% calculates the gradient at any point on the RBF surface

if (nargout > 1)
    
    g = zeros(d,1);
    
    phiprime_vector = lambda.*(phiprime(norm_value,rbftype,gamma)./(norm_value+realmin));
        
    for j = 1:d
        
        tempsum = phiprime_vector'*(xtemp(:,j)-D(:,j)+realmin);
        
        switch polytail
        case 'linear'
            tempsum = tempsum + ctail(1+j);
        case 'quadratic'
            tempsum = tempsum + ctail(1+j);
            pos = 1 + d + j;
            for k = 1:(j-1)
                tempsum = tempsum + ctail(pos)*x(k);
                pos = pos + (d-k);
            end
            tempsum = tempsum + 2*ctail(pos)*x(j);
            for k = (j+1):d
                pos = pos + 1;
                tempsum = tempsum + ctail(pos)*x(k);
            end
        end
                
        g(j) = tempsum;
        
    end
    
end



% computes initial matrices for RBF interpolation

function [PHI, phi0, P, pdeg, pdim] = ...
    ComputeInitialRBFMatrices(maxsize, numinitpts, rbftype, gamma, polytail)

global d D;

% determine pairwise distance between points
PairwiseDistance = zeros(numinitpts);
for i = 1:(numinitpts-1)
    xtemp = ones(numinitpts-i, 1)*D(i, :);
    PairwiseDistance(i, i+1:numinitpts) = sqrt(sum((xtemp-D(i+1:numinitpts, :)).^2, 2))';
end

% initial PHI matrix
PHI = zeros(maxsize);
phi0 = phi(0, rbftype, gamma);
for i = 1:numinitpts
    PHI(i, i) = phi0;
    PHI(i, i+1:numinitpts) = phi(PairwiseDistance(i, i+1:numinitpts), rbftype, gamma);
    PHI(i+1:numinitpts, i) = PHI(i, i+1:numinitpts)';
end

% initial P matrix
switch polytail
case 'none'
    pdeg = -1;
    pdim = 0;
    P = [];
case 'constant'
    pdeg = 0;
    pdim = 1;
    P = ones(maxsize, 1);
case 'linear'
    pdeg = 1;
    pdim = d + 1;
    P = [ones(maxsize, 1), D];
case 'quadratic'
    pdeg = 2;
    pdim = ((d+1)*(d+2))/2;
    P = [ones(maxsize, 1), D, zeros(maxsize, (d*(d+1))/2)];
    columnpos = d+1;
    for i = 1:d
        for j = i:d
            columnpos = columnpos + 1;
            P(1:numinitpts, columnpos) = D(1:numinitpts, i).*D(1:numinitpts, j);
        end
    end
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end



% sets LMSRBF internal algorithm parameters

function SetLocalStochRBFAlgorithmParams

global d minxrange maxsigma_stdev minsigma_stdev sigma_stdev numcand ...
    weightfactors numweightfactors;

maxsigma_stdev = 0.1*minxrange;             % maximum step size
minsigma_stdev = 0.0005*minxrange;          % minimum step size
sigma_stdev = maxsigma_stdev;               % initial step size
numcand = min(1000*d, 20000);               % number of normal random candidates points in each iteration
weightfactors = [0.95];                     % weights for the RBF criterion (another option: weightfactors = [0.5, 0.95])
numweightfactors = length(weightfactors);   % number of weight factors



% initializes LMSRBF control parameters

function InitializeLocalStochRBFControlParams

global d minxrange successthreshold failtolerance successctr failctr ...
    innerballradius outerballradius numreqininnerball numreqinouterball localminflag;

successthreshold = 3;                   % threshold for the number of consecutive successes
failtolerance = min(max(d,5),30);       % threshold for the number of consecutive failures
successctr = 0;                         % number of consecutive successful iterations
failctr = 0;                            % number of consecutive unsuccessful iterations
innerballradius = 0.0025*minxrange;     % radius of inner trust region ball at xbest
outerballradius = 2*innerballradius;    % radius of outer trust region ball at xbest 
numreqininnerball = d+1;                % number of evaluated points required in inner ball centered at xbest
numreqinouterball = 2*numreqininnerball;% number of evaluated points required in outer ball centered at xbest
localminflag = 0;                       % indicates whether or not xbest is at a local minimum



% replaces function values by the median of all available function
% values if they exceed the median
% transform is either 'T' (with transformation) or 'N' (no transformation)
% Ftransform is the vector of possibly transformed function values

function Ftransform = TransformFunctionValues(transform)

global n F;

Ftransform = F(1:n);
if (strcmp(transform,'T') == 1)
    medianF = median(F(1:n));
    for i = 1:n
        if (Ftransform(i) > medianF)
            Ftransform(i) = medianF;
        end
    end
end



% fits an RBF model using all available data points
% lambda is the vector of coefficients for the RBF part
% ctail is the vector of coefficients for the polynomial tail

function [lambda, ctail] = FitRBF

global polytail pdim n PHI P Ftransform;

if (strcmp(polytail,'none') == 1)
    coeff = PHI(1:n,1:n)\Ftransform;
else
    Ptemp = P(1:n,:);
    coeff = [PHI(1:n,1:n), Ptemp; Ptemp', zeros(pdim)]\[Ftransform; zeros(pdim,1)];
end
lambda = coeff(1:n);
ctail = coeff(n+1:n+pdim);



% selects LMSRBF iterate

function [xselected, normval] = SelectLocalStochRBFPoint(iterctr)

global d xmin xmax tolerance xbest sigma_stdev numcand weightfactors numweightfactors;

% set the weight of the RBF criterion
weightindex = mod(iterctr, numweightfactors);
if (weightindex == 0)
    weightindex = numweightfactors;
end
valueweight = weightfactors(weightindex);
% fprintf('RBF weight = %f\n',valueweight);

% generate candidate points
fprintf('step size = %f\n',sigma_stdev);
CandPoints = GenerateBoundedNormalRandPoints(xbest, sigma_stdev*ones(1,d), xmin, xmax, numcand);

% compute the total weighted score
% [CandTotalValues, NormValues, ScaledCandValues, ScaledCandMinDist, CandValues, CandMinDist] = ...
%     ComputeStochRBFWeightedScore(CandPoints, valueweight, tolerance);
[CandTotalValues, NormValues] = ComputeStochRBFWeightedScore(CandPoints, valueweight, tolerance);

% select the best candidate point for function evaluation
% [MinCandTotalValues,selindex] = min(CandTotalValues);
[~,selindex] = min(CandTotalValues);
xselected = CandPoints(selindex,:);
normval = NormValues(:,selindex)';

% % info on selected point
% fprintf('ScaledCandValue for selected point: %f\n',ScaledCandValues(selindex));
% fprintf('ScaledCandMinDist for selected point: %f\n',ScaledCandMinDist(selindex));
% fprintf('CandTotalValue for selected point: %f\n',CandTotalValues(selindex));



% generates normal random points (covariance matrix is assumed to be diagonal)
% xcenter is the mean, stdev is the vector of standard deviations
% points generated outside the domain [xmin, xmax] will be put on the boundary 

function RandPoints = GenerateBoundedNormalRandPoints(xcenter, stdev, xmin, xmax, numrand)

dim = length(xmin);
RandPoints = zeros(numrand, dim);
for i = 1:numrand
    RandPoints(i,:) = max( xmin, min(xcenter + stdev.*randn(1,dim), xmax) );
end



% computes the total weighted score for each candidate point
%
% Inputs: 
% (1) CandPoints is a matrix whose rows are the candidate points
% (2) valueweight is the weight for the RBF criterion
%
% Outputs:
% (1) CandTotalValues is a vector of total weighted scores for the 
% candidate points
% (2) NormValues is a matrix of distances between the previously evaluated
% points and the candidate points
% (3) ScaledCandValues is a vector of scores for the RBF criterion
% (4) ScaledCandMinDist is a vector of scores for the distance criterion

function [CandTotalValues, NormValues, ScaledCandValues, ScaledCandMinDist, CandValues, CandMinDist] = ...
    ComputeStochRBFWeightedScore(CandPoints, valueweight, distreq)

% estimate the function value of the candidate points
[CandValues, NormValues] = ComputeRBF(CandPoints);
MinCandValues = min(CandValues);
MaxCandValues = max(CandValues);
ScaledCandValues = (CandValues-MinCandValues)/(MaxCandValues-MinCandValues);

% determine the minimum distance of candidate points from previously evaluated points
CandMinDist = (min(NormValues,[],1))';
MaxCandMinDist = max(CandMinDist);
MinCandMinDist = min(CandMinDist);
ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);

% compute the total weighted score
mindistweight = 1 - valueweight;
CandTotalValues = valueweight*ScaledCandValues + mindistweight*ScaledCandMinDist;

% remove candidate points that are too close to previously evaluated points
numcand = size(CandPoints, 1);
numdiscard = sum(CandMinDist < distreq);
if (numdiscard < numcand)
    CandTotalValues(CandMinDist < distreq) = inf;
end
%fprintf('Number of discarded points (too close to a previously evaluated point): %d\n',numdiscard);



% Compute the value of the RBF surface for several points at once
%
% syntax:
% RBFVALUE = ComputeRBF(Ymatrix)

function [RBFVALUE, NORMVALUE, U_Y] = ComputeRBF(Ymatrix)

global d n D rbftype gamma polytail lambda ctail;

numpoints = size(Ymatrix,1);

% Ymatrix is an numpoints x d matrix
% size of D and lambda vary although D always has d columns

NORMVALUE = zeros(n,numpoints);
for k = 1:numpoints
    ytemp = ones(n,1)*Ymatrix(k,:);
    NORMVALUE(:,k) = sqrt(sum((ytemp-D(1:n,:)).^2,2));
end
U_Y = phi(NORMVALUE,rbftype,gamma);

switch polytail
case 'none'
    PolyPart = zeros(numpoints,1);
case 'constant'
    PolyPart = ctail*ones(numpoints,1);
case 'linear'
    PolyPart = [ones(numpoints,1),Ymatrix]*ctail;
case 'quadratic'
    temp = [ones(numpoints,1),Ymatrix,zeros(numpoints,(d*(d+1))/2)];
    columnpos = d+1;
    for i = 1:d
        for j = i:d
            columnpos = columnpos + 1;
            temp(:,columnpos) = Ymatrix(:,i).*Ymatrix(:,j);
        end
    end
    PolyPart = temp*ctail;
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end

RBFVALUE = (U_Y')*lambda + PolyPart;



% updates the control parameters for the local stochastic RBF algorithm

function UpdateLocalStochRBFControlParams

global D n xbest Fbest Fbest_old maxsigma_stdev minsigma_stdev sigma_stdev ...
    successthreshold failtolerance successctr failctr innerballradius outerballradius ...
    numreqininnerball numreqinouterball localminflag;

if (Fbest < Fbest_old)
    successctr = successctr + 1;    
    failctr = 0;
%     disp('Improvement in function value!');
else
    successctr = 0;    
    failctr = failctr + 1;
end

if (successctr >= successthreshold)
    
    % increase step size
    sigma_stdev = min(2*sigma_stdev, maxsigma_stdev);
    fprintf('Doubled the step size! (new step size = %f)\n',sigma_stdev);
    
    % reset successctr
    successctr = 0;
    
end

if (failctr >= failtolerance)
        
    % reset failctr
    failctr = 0;
        
    % reduce step size
    sigma_stdev = max(sigma_stdev/2, minsigma_stdev);
    fprintf('Reduced step size by a half! (new step size = %f)\n',sigma_stdev);

end


% check if algorithm is in a local minimum
xbestdist = ComputeDistanceMatrix(xbest, D(1:n,:));
numevalptsininnerball = sum(xbestdist<=innerballradius);    % number of points in a inner ball centered at xbest
numevalptsinouterball = sum(xbestdist<=outerballradius);    % number of points in a outer ball centered at xbest
% disp('-------------------------------------------------------------------------');
% fprintf('Number of points in inner ball centered at the current best solution: %d\n',numevalptsininnerball);
% fprintf('Number of points in outer ball centered at the current best solution: %d\n',numevalptsinouterball);
% disp('-------------------------------------------------------------------------');
if ((numevalptsininnerball >= numreqininnerball) && (numevalptsinouterball >= numreqinouterball))
    % [validflag, GoodPoints, NFP] = CheckModel(xbest, innerballradius);
    [validflag, ~, ~] = CheckModel(xbest, innerballradius);
    if (validflag == 1)
        localminflag = 1;
        disp('Current best solution is near a local minimum!');
        disp('Restarting LMSRBF using a new symmetric Latin hypercube design (SLHD).');
    end    
end



% computes the distances between two sets of points
% X and Y are matrices whose rows are points
% Distances is the matrix of distances between points in X and points in Y

function Distances = ComputeDistanceMatrix(X, Y)

numxrows = size(X, 1);
numyrows = size(Y, 1);
Distances = zeros(numxrows, numyrows);
for i = 1:numxrows
    Xtemp = ones(numyrows, 1)*X(i,:);
    Distances(i,:) = sqrt(sum((Xtemp-Y).^2, 2))';
end



% checks whether the RBF model is valid at a given point

function [validflag, GoodPoints, NFP] = CheckModel(xcurrent, TRradius)

% global d xmin xmax n D;
global d n D;

% additional parameters
C1 = 2;
theta1 = 0.001;

%%%%% validity step %%%%%

% initialize the intermediate polynomials
NFP = [zeros(d,1), eye(d)];
AssignIndicator = zeros(d,1);

% collect the good points
GoodPoints = zeros(d+1, d);

% select the current iterate
GoodPoints(1,:) = xcurrent;

% update the NFPs
NFP(:,1) = NFP(:,1) - NFP*[1,xcurrent]';

% collect points inside the trust region
xcurrentInfDist = ComputeInfNormDistanceMatrix(xcurrent, D(1:n,:));
TRpointindicator = xcurrentInfDist<=C1*TRradius;
numTRpoints = sum(TRpointindicator);
TRpoints = D(TRpointindicator,:);
TRinfdist = xcurrentInfDist(TRpointindicator);
% [minTRinfdist,index] = min(TRinfdist);  % index of xcurrent in TRpoints
[~,index] = min(TRinfdist);  % index of xcurrent in TRpoints

% initialize indicator of selected points
SelectedIndicator = zeros(numTRpoints,1);
SelectedIndicator(index) = 1;   % xcurrent has been selected

% compute pivots
numgoodpoints = 1;
for i = 1:d
    pivots = zeros(numTRpoints,1);
    for j = 1:numTRpoints
        if (SelectedIndicator(j) == 0)
            pivots(j) = NFP(i,:)*[1,TRpoints(j,:)]';
        end
    end
    [maxpivot, selindex] = max(abs(pivots));
    if (maxpivot > theta1)
        
        % update set of good points
        numgoodpoints = numgoodpoints + 1;
        GoodPoints(numgoodpoints,:) = TRpoints(selindex,:);
        SelectedIndicator(selindex) = 1;
        AssignIndicator(i) = 1;
        
        % normalize current polynomial
        NFP(i,:) = NFP(i,:)/pivots(selindex);
        
        % update polynomials
        for j = 1:d
            if (j ~= i)
                NFP(j,:) = NFP(j,:) - (NFP(j,:)*[1,TRpoints(selindex,:)]')*NFP(i,:);
            end
        end
        
    end
end

if (numgoodpoints < d+1)
    validflag = 0;
else
    validflag = 1;
end



% computes the distances between two sets of points
% X and Y are matrices whose rows are points
% Distances is the matrix of distances between points in X and points in Y

function Distances = ComputeInfNormDistanceMatrix(X, Y)

numxrows = size(X, 1);
numyrows = size(Y, 1);
Distances = zeros(numxrows, numyrows);
for i = 1:numxrows
    Xtemp = ones(numyrows, 1)*X(i,:);
    Distances(i,:) = (max(abs(Xtemp-Y),[],2))';
end



% updates the RBF matrices PHI and P after a single function evaluation
% xselected is the newly evaluated point
% normval is the vector of distances between xselected and all previously evaluated points

function UpdateRBFMatrices(xselected, normval)

global d rbftype gamma polytail n PHI phi0 P;

n_old = n - 1;

% compute new_phi
new_phi = phi(normval, rbftype, gamma);

% update PHI and P
PHI(n,1:n_old) = new_phi;
PHI(1:n_old,n) = new_phi';
PHI(n,n) = phi0;

switch polytail
    case 'linear'
        P(n,2:d+1) = xselected;
    case 'quadratic'
        P(n,2:d+1) = xselected;
        columnpos = d+1;
        for i = 1:d
            for j = i:d
                columnpos = columnpos + 1;
                P(n,columnpos) = xselected(i)*xselected(j);
            end
        end
end



% calculates the value of the radial function

function output = phi(r,type,gamma)

switch type
case 'linear'
    output = r;
case 'cubic'
    output = r.^3;
case 'thinplate'
    if (r >= 0)
        output = (r.^2).*log(r+realmin);
    else
        output = zeros(size(r));
    end
case 'multiquadric'
    output = sqrt(r.^2 + (gamma*ones(size(r))).^2);
case 'Gaussian'
    output = exp(-(r./(gamma*ones(size(r)))).^2);
otherwise
    disp('Error: Unknown type.');
    return;
end



% calculates the derivative of the radial function

function output = phiprime(r,type,gamma)

switch type
case 'linear'
    output = ones(size(r));
case 'cubic'
    output = 3*(r.^2);
case 'thinplate'
    if (r >= 0)
        output = r + 2*(r.*log(r+realmin));
    else
        output = zeros(size(r));
    end
case 'multiquadric'
    output = r./sqrt(r.^2 + (gamma*ones(size(r))).^2);
case 'Gaussian'
    output = -((2*r)./((gamma*ones(size(r))).^2)).*exp(-(r./(gamma*ones(size(r)))).^2);
otherwise
    disp('Error: Unknown type.');
    return;
end



% gets the final output values for the algorithm
%
% Outputs:
% (1) solution is the best solution obtained by the algorithm
% (2) value is the best function value obtained (i.e. function value of the
% best solution)
% (3) numevals is the total number of function evaluations
% (4) proctime is the total processing time
% (5) totalfevaltime is the total time spent on all function evaluations
% (6) success indicates whether the algorithm reached the specified
% accuracy in the solution (success=1) or not (success=0) 

function [solution, value, numevals, proctime, totalfevaltime, success] = GetOutputValues

global targetvalue fevaltime n time0 xbest Fbest;

solution = xbest;
value = Fbest;
numevals = n;
proctime = etime(clock, time0);
totalfevaltime = sum(fevaltime(1:n));
if (Fbest <= targetvalue)
    success = 1;
else
    success = 0;
end