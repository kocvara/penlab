% test file for penlabtest.m

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 29 August 2019

function penm = test_Lyapunov()

A=[0 1 0;0 0 1;-1 -2 -3];
penm = lyapu(A);
