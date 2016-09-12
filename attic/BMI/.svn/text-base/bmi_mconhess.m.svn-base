function [Akddx, userdata] = bmi_mconhess(x,Y,k,i,j,userdata)
% Compute derivatives: d/dx_i A_k(x) based on the data from sdpdata
% in this context it is in fact -F_i of the specific block

  Akddx=[];
  if (k<=0 || k>userdata.Na)
    return;
  end
  if (i<=0 || i>userdata.Nx)
    return;
  end
  if (j<=0 || j>userdata.Nx)
    return;
  end

  if i==length(x), return; end
  indK=userdata.indK{k};
  li = indK(i,1);
  if li>0
      for jj=1:li
          j1 = indK(1,jj+1);
          if j1==j
              Akddx = - (userdata.K{k,i,j1});
          end
      end
  end


