clear
ds = 1; 
load('avgcg.mat');
c0 =reshape(avgcg,[3 3 3 3]); 

c66 =MS_cijkl2cij(c0);
s66 = pinv(c66);
s0 = MS_cij2cijkl(s66); 



nn=[1024 1024 4];
N1 = nn(1);
N2=nn(2);
N3= nn(3);
for kxx= 1: ds: N1
    for kyy=1: ds :N2
        for kzz=1: ds :N3
             
                    
                    
     % compute gamma
[   Gamma_fft  xk] = compute_gamma_fft( kxx, kyy, kzz, c0,s0, nn);
 
       %if norm(xk)==0  %0 frequency case
            
        %    Gamma_fft = -s0;
       %end
 for i = 1:3
for  j = 1:3
 for k=1:3
for l=1:3            
 Gfft(i,j,k,l,kxx, kyy,kzz) = Gamma_fft(i,j,k,l);
end
end
end
end
            

        end
    end
end

disp('writing')
for ii =1:3
for jj=1:3
for kk=1:3
for ll= 1:3
G= Gfft(ii,jj,kk,ll,:,:,:);
G=squeeze(G);
matn = sprintf('Gfft%d%d%d%d.mat', ii,jj,kk,ll)
save(matn,'G')
G=[];
end
end
end
end

