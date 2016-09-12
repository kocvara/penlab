function [dg, userdata]=bmi_congrad(x,Y,userdata)
% returns all inequalities at once, expect g(x)<=0
% vector Ng x Nx

  if (userdata.Ng>0)
    dg=userdata.B';
  else
    dg=[];
  end

