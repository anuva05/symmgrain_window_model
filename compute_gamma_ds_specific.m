function Gfft = compute_gamma_ds_specific(N1, N2, N3, c0, s0, nn,ei,ej, ii,jj, ds )

for kxx= 1: ds: N1
    for kyy=1: ds :N2
        for kzz=1: ds :N3
             
                    
                    
     % compute gamma
[   Gamma_fft  xk] = compute_gamma_fft( kxx, kyy, kzz, c0,s0, nn);
 
  
       %if norm(xk)==0  %0 frequency case
            
        %    Gamma_fft = -s0;
       %end
            
  i=ei; %component of strain being computed
  j = ej;
  k=ii;
  l=jj;
 Gfft((kxx+(ds -1))/ds, (kyy+(ds -1))/ds, (kzz+(ds -1))/ds ) = Gamma_fft(i,j,k,l); %only store 1 comp
             

        end
    end
end
