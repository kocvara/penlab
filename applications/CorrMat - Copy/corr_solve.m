% Various instances of the nearest correlation matrix problem

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 19 August 2016

kappa=10;
%% nearest correlation matrix
penm = corr_define;

%% nearest correlation matrix with element-wise constraints
penm = corr_define_bound;

%% nearest correlation matrix with constrained condition number
%penm = corr_define_cond(kappa);

%% nearest correlation matrix with constrained condition number and elements
%penm = corr_define_cond_bound(kappa);

problem = penlab(penm);
problem.solve();

%% no constraints on condition number
X = problem.Y{1}

%% constraints on condition number
%X = problem.Y{1}*problem.x
