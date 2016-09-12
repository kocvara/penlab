function [Akdx, userdata] = bmi_mcongrad(x,Y,k,i,userdata)
% Compute derivatives: d/dx_i A_k(x) based on the data from sdpdata
% in this context it is in fact -F_i of the specific block

  Akdx=[];
  if (k<=0 || k>userdata.Na)
    return;
  end
  if (i<=0 || i>userdata.Nx)
    return;
  end

  Akdx=-userdata.A{k,i+1};
  
  indK=userdata.indK{k};

  if i==length(x), return; end
  
  li = indK(i,1);
  if li>0
      for jj=1:li
          j = indK(1,jj+1);
          Akdx = Akdx - x(j).*userdata.K{k,i,j};
      end
  end
  
  for ii = 1:i-1
      li = indK(ii,1);
      if li>0
          for jj=1:li
              j = indK(1,jj+1);
              if j==i
                  Akdx = Akdx - x(ii).*userdata.K{k,ii,j};
              end
          end
      end
  end


