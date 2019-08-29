function penm = test_TTOB()
%
% test file for penlabtest.m

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 29 August 2019

warning('off','MATLAB:nearlySingularMatrix')
par = kobum('GEO/t3x3.geo');
penm = tto_define(par);


