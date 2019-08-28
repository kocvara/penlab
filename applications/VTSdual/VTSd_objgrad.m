function [df, userdata]=VTSd_objgrad(x,Y,userdata)

  nelem = userdata.nelem;
  nnod = userdata.nnod;
  rhs = userdata.rhs;
  ubound = userdata.ubound;
  V = userdata.V;
  
  alpha = x(nnod+1);
  
  df = [-rhs; V; ubound.*ones(nelem,1)];

