%function  solve_fmo
%
% Solve the free material optimization problem by PENLAB
%

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013


par = readdata;
problem = fmo_define(par);
penm = penlab(problem);
penm.opts.mpenalty_update=.5;
%penm.opts.mmlt_update=.3;
%penm.opts.inner_stop_limit=.0001;
%penm.opts.max_outer_iter=7;
%penm.opts.outlev=3;
solve(penm);

% figure plot
nelx=par.nelx;
nely=par.nely;
nelem = par.nelem;
for i=1:nelem
    trE(i) = trace(penm.Y{i});
end

aa = reshape(trE(1:nelem),nelx,nely);
imagesc(aa'); axis equal; axis tight; axis off; set(gca, 'YDir', 'normal');

%end

