% *Linear semidefinite programming with SDPA input*
% 
% Assume that a linear SDP problem is stored in an SDPA input file, for
% instance |arch0.dat-s| stored in directory |datafiles|. All the user has
% to do to solve the problem is to call the following sequence of commands:

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 27 Nov 2013

sdpdata=readsdpa('../../datafiles/control1.dat-s');
penm=sdp_define(sdpdata)
prob=penlab(penm);
prob.solve();