function  penm = test_VTSdual()

% test file for penlabtest.m

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 29 August 2019
%
load stiffmat2   % stiffmat1 ... stiffmat5 

par.A = A; par.IA = IA; par.nelem=nelem; par.nnod=nnod;

if nelem==16
    nx=4;ny=4;        %stiffmat1
elseif nelem==200
    nx=20;ny=10;      %stiffmat2
elseif nelem==420
    nx=30;ny=14;      %stiffmat3
elseif nelem==800
    nx=40;ny=20;      %stiffmat4
elseif nelem==1800
    nx=60;ny=30;      %stiffmat5
elseif nelem==3200
    nx=80;ny=40;      %stiffmat6
elseif nelem==5000
    nx=100;ny=50;      %stiffmat7
elseif nelem==7200
    nx=120;ny=60;      %stiffmat8    
else
    display('error in nelem');
end

par.nx = nx; par.ny = ny;

par.nloads = 1;

RHS = circshift(RHS,1);
par.ff{1} = RHS;

penm = VTSdual_lu_define(par);


