function [ddf, userdata]=sVTSd_congrad(x,Y,k,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  ubound = userdata.ubound;
  A = userdata.A;
  
  ddf = sparse(nelem+1+nnod,nelem+1+nnod);
  ddf(1:nnod,1:nnod) = A{1,k};


