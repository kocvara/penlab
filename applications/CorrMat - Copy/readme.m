% We solve the nearest correlation matrix problem. The problem may include
% constraints on matrix elements and on the condition number of the
% resulting matrix; see Example 7.1 from the PENLAB paper (directory
% tex/penlab_paper). On input is a given matrix H, not necessarily a
% correlation matrix and, if applicable, kappa, the required condition
% number of the computed correlation matrix.

% Choose (by commenting/uncommenting) the required problem in |corr_solve|:
% nearest correlation matrix
penm = corr_define;
% nearest correlation matrix with element-wise constraints
penm = corr_define_bound;
% nearest correlation matrix with constrained condition number
penm = corr_define_cond(kappa);
% nearest correlation matrix with constrained condition number and elements
penm = corr_define_cond_bound(kappa);

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013-2016 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 13 Sep 2016