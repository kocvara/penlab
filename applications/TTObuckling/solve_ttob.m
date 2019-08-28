function solve_ttob( iname )
%
% Solve the truss topology optimization problem with buckling constraint by
% PENLAB
%
% Example: >> solve_ttob('GEO/t3x3.geo')

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013
warning('off','MATLAB:nearlySingularMatrix')
par = kobum(iname);
%par.f(8)=0;par.f(7)=1;
problem = tto_define(par);
penm = penlab(problem);
solve(penm);
pic(par,penm.x);

end

