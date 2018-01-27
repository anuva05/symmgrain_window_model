function[ Gamma_fft  xk] = compute_gamma_fft( kxx, kyy, kzz, c0,s0, nn)
 
 
delt = [1 1 1];
npts1 = nn(1);
npts2 = nn(2);
npts3 = nn(3);
 
      %%%%%% compute greeens %%%%%%
       
       if(kxx<=npts1/2)
          kx=kxx-1;
      end
      
      if(kxx>npts1/2)
          kx=kxx-npts1-1;
      end
      
 
      if(kyy<=npts2/2)
          ky=kyy-1;
      end
      if(kyy>npts2/2)
          ky=kyy-npts2-1;
      end
      
 
      if(kzz<=npts3/2)
          kz=kzz-1;
      end
      if(kzz>npts3/2)
          kz=kzz-npts3-1;
      end
    
 
      xk(1)=kx/(delt(1)*nn(1));
      xk(2)=ky/(delt(2)*nn(2));
 
      xk(3)=kz/(delt(3)*nn(3));
 
 
      xknorm=sqrt(xk(1)^2+xk(2)^2+xk(3)^2);
 
      if (xknorm~=0) 
      for  i=1:3
      xk2(i)=xk(i)/(xknorm*xknorm*2.*pi);
      xk(i)=xk(i)/xknorm;
      end
      end
      
 
 
      if(kxx==npts1/2+1)|(kyy==npts2/2+1)| (kzz==npts3/2+1) 
         Gamma_fft=-s0;
      else
      for i=1:3
          for k=1:3
                    a(i,k)=0;
              for j=1:3
                  for l=1:3
  
      a(i,k)=a(i,k)+ c0(i,j,k,l)*xk(j)*xk(l);
                  end
              end
          end
      end
  
      
 
 
      a = pinv(a);
      
 
 
% the inverse of A replaces A
  for p=1:3
          for q=1:3
 
              for i=1:3
                  for j=1:3
      Gamma_fft(p,q,i,j)=-a(p,i)*xk(q)*xk(j);
                end
              end
          end
  end
       
                end %if loop ends 
