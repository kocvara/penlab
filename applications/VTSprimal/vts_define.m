function [penm] = tto_define(par)
% define penm structure for the VTS primal problem

% This file is a part of PENLAB package distributed under GPLv3 license
% Copyright (c) 2013 by  J. Fiala, M. Kocvara, M. Stingl
% Last Modified: 30 Aug 2019

nx = par.nx; ny = par.ny;
A = par.A; IA = par.IA; nelem=par.nelem; nnod=par.nnod;
nloads = par.nloads;
ny2 = ny/2; nelem2 = nelem/2; nnod2 = 2*nx*ny2; nnodi = 2*nx;
n1 = nnod ; m = nelem;
volume = nelem/3;

lb = 1;
for ie=1:nelem
    len = IA(ie);
    %clear ss;
    ss = sparse(A(lb:lb+len-1),A(lb+len:lb+2*len-1),A(lb+2*len:lb+3*len-1),nnod,nnod);
    hel = triu(ss) + tril(ss)' - diag(diag(ss));
    ss = hel +hel' - diag(diag(hel)) ;
    AA{1,ie+1} = [0 sparse(1,n1); sparse(n1,1) ss];
    lb = lb + 3*len;
end


AA{1,m+2} = sparse(n1+1,n1+1);
AA{1,m+2}(1,1) = 1;

for k=1:nloads
    ff = par.ff{k};
    AA{k,1} = -[0 -ff'; -ff sparse(n1,n1)];
end

sdpdata.Nx = m+1;% ..... number of primal variables
sdpdata.Na = nloads;% ..... number of linear matrix inequalities
sdpdata.Ng = 1;% ..... number of linear inequalitites
sdpdata.c = [zeros(m,1);1];% ... dim (Nx,1), coefficients of the linear objective function
sdpdata.NaDims = [n1+1];% . vector of sizes of matrix constraints (diagonal blocks)
sdpdata.A = AA;
sdpdata.Adep = 1:m+1; ... rather create here vvv !!!
    clear A;

penm = [];


% keep the whole structure
penm.userdata=sdpdata;

penm.Nx=sdpdata.Nx;
%penm.lbx=[zeros(m,1);-Inf];
penm.lbxbar=[zeros(m,1);-Inf];
%penm.ubx=[1.7.*ones(m,1);Inf];
penm.ubx=[1.0.*ones(m,1);Inf];

penm.NgLIN=sdpdata.Ng;  %length(sdpdata.d);
penm.ubg = [volume];

penm.NALIN=sdpdata.Na;
% let's make it negative semidefinite
penm.ubA=zeros(sdpdata.Na,1);

penm.objfun = @vts_objfun;
penm.objgrad = @vts_objgrad;
penm.objhess = @vts_objhess;

penm.confun = @vts_confun;
penm.congrad = @vts_congrad;
%penm.conhess = @sdp_conhess;  not needed because all linear

penm.mconfun = @vts_mconfun;
penm.mcongrad = @vts_mcongrad;
% hessian not needed as linear

end

function [g,userdata] = vts_confun(x,Y,userdata)

g=sum(x(1:end-1));

end

function [f,userdata] = vts_objfun(x,Y,userdata)

f = userdata.c'*x;

end

function [ddf, userdata] = vts_objhess(x,Y,userdata)

ddf = [];

end

function [Akx, userdata] = vts_mconfun(x,Y,k,userdata)
% evaluate A_k(x) based on sdpdata, k denotes a block number
% sdpdata is a structure as obtained from readsdpa.m
% note that we aim for A_k(x)<=0 thus A_k(x) is 'reversed' than usual:
%   A_k(x) = F_0 - sum x_i*F_i


Akx=[];
if (k<=0 || k>userdata.Na)
    return;
end

Akx=userdata.A{k,1};
for i=1:userdata.Nx
    % won't work if F{} is [] ... <== dims must match
    Akx = Akx - x(i).*userdata.A{1,i+1};
    % or - ?
end

end

function [Akdx, userdata] = vts_mcongrad(x,Y,k,i,userdata)

Akdx=[];
if (k<=0 || k>userdata.Na)
    return;
end
if (i<=0 || i>userdata.Nx)
    return;
end

Akdx=-userdata.A{1,i+1};

end

function [df, userdata]=vts_objgrad(x,Y,userdata)

df = userdata.c;

end

function [dg, userdata]=vts_congrad(x,Y,userdata)

dg = [ones(length(x)-1,1);0];

end




