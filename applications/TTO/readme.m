% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013-2016 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 30 Aug 2019

% We solve the truss topology optimization problem, primal formulation,
% written as linear SDP.

% min_x sum(x)
% subject to
% [compl ff'; ff zero(n,n)] + sum_i [0 zero(1,n); zero(n,1) A_i)] \succeq 0;

% where compl \in \R, ff \in \RR^n, A_i \in \RR^{n\times n} sym, psd, i=1...m,
% are given and x \in \RR^m is the variable

% Use

solve_tto('GEO/t3x3.geo');

% to solve a sample problem. Here the initial ground structure t3x3.geo 
% can be replaced by any ground structure from directory GEO

