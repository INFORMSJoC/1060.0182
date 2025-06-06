% Interface for running the MLMSRBF algorithm
%
% MLMSRBF.p implements the Multistart Local Metric Stochastic RBF
% (MLMSRBF) method for finding the global minimum of a function defined
% over a hypercube, i.e., the only constraints are bound constraints on all
% the decision variables. The method is described in the paper:
%
% R.G. Regis, C.A. Shoemaker. A stochastic radial basis function method 
% for the global optimization of expensive functions. INFORMS Journal on
% Computing, Vol. 19, No. 4, pp. 497-509, 2007.
%
% Written by: Rommel G. Regis
% Last updated: 07/03/07, 06/30/10, 07/05/10
%
% Syntax:
% outputinfo = RunMLMSRBF(probinfo, algparams)
%
%
% Input:
%
%   probinfo.funcname: name of function represented as a string
%
%   probinfo.lb: row vector of lower bounds (does not allow -inf)
%
%   probinfo.ub: row vector of upper bounds (does not allow inf)
%
%   probinfo.lb and probinfo.ub define the bound constraints for the
%   problem. The method requires that the search space is a hypercube so
%   this means that the difference (probinfo.ub - probinfo.lb) must be a
%   constant vector whose constant entry is strictly positive. Ideally, we
%   recommend that the problem be transformed so that probinfo.lb is the
%   zero vector and probinfo.ub is the all ones vector.
%   
%   (OPTIONAL) probinfo.userinitpoints: matrix of initial points provided
%   by the user (The initial points are the rows of this matrix.) (default
%   value = [])
%
%   (OPTIONAL) algparams.rbftype: type of RBF interpolant (Set to 'cubic'
%   or 'thinplate') (default value = 'thinplate')
%                             
%   algparams.maxeval: maximum number of function evaluations 
%   (The algorithm terminates when the total number of calls to the
%   objective function exceeds maxeval.)
%
%   (OPTIONAL) algparams.numdesignpts: number of space-filling experimental
%   design points for every run of LMSRBF (must be >= dim+1, where dim is
%   the number of decision variables) (default value = 2*(dim + 1))
%
%   (OPTIONAL) algparams.transform: indicates whether the high function 
%   values will be replaced by the median of all available function values
%   before fitting the RBF interpolation model (Set to 'T' to replace high
%   function values or 'N' to use original function values.) 
%   (default value = 'T')
%
%
% Output:
%
%   outputinfo.bestsolution: best solution vector found
%
%   outputinfo.bestvalue: function value of best solution vector (i.e., the
%   best function value found
%
%   outputinfo.numevals: actual number of calls to the objective function
%
%   outputinfo.history: history matrix that summarizes the entire history 
%   and trajectory of the algorithm from the beginning until termination.
%   It is a matrix of size numevals x (dim+4), where numevals is the number
%   of calls to the objective function and dim is the dimension of the
%   problem (or number of decision variables). Each row of the history
%   matrix corresponds to one function evaluation. In particular, for the
%   nth function evaluation, we have:
%
%   history(n,1) is the LMSRBF run number
%   history(n,2) is the function value of the nth function evaluation point
%   history(n,3:dim+2) is the nth function evaluation point
%   history(n,dim+3) is the time (in seconds) to execute the nth function
%   evaluation
%   history(n,dim+4) is time elapsed (in seconds) up to the nth function
%   evaluation excluding the time spent on all function evaluations so far
%   (i.e., excluding the total time for the first n function evaluations).

function outputinfo = RunMLMSRBF(probinfo, algparams)

outputinfo = MLMSRBF(probinfo, algparams);