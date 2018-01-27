%fiber in  + out model error
%single fiber modeling and replacing with approximation in fortran
%iterations
%clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maxNumCompThreads
%mat_file=('sq128_symmgrain.mat');
%load(mat_file);
%infile = 'C:\Users\anuvak\Documents\fortran_data_process\SG10131.out';
%infile;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% get s0,c0,nn,s11 from files written out
s = dlmread('s0.out');
s1 = s(:);
s1 =s1(1:end-4);
s0 = reshape(s1,[3 3 3 3]); 
nn = dlmread('nn.out');
c = dlmread('c0initial.out');
c1 = c(:);
c1 = c1(1:end-4);
c0= reshape(c1,[3 3 3 3]);


 save init.mat nn c0 s0
% NOT SURE IF S0 AND C0 SHOULD BE FOR EACH POINT OR AVERAGE TO COMPUTE GREENS?? 

% read al sg
for kz=1:nn(3)
   for i=1:3
      for j=1:3
               sfile=sprintf('SG%d%d%d.out', (kz+100),i,j);   
               sp = dlmread(sfile);
              %  norm(sp) checking if norm is nan
               splane = reshape(sp,nn(1),nn(2));
               sgfull(i,j,:,:,kz) = splane;
      end 
     end
end





%%%%%%%%% set %%%%%%%%%%%%%%%%%%%%%%%%%%

ds_in =2 ;%downsample factor
type_win=1; %window type.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaE =zeros(3,3,nn(1),nn(2),nn(3));
deltaE2 =zeros(3,3,nn(1),nn(2),nn(3));
for ei=1:3
  for ej=1:3


        for ii=1:3
        for jj=1: 3
         sc = squeeze(sgfull(ii,jj,:,:,:));

        [total_interp_out, recon_interp_out]  = compress_update_pipeline2d(c0,s0,nn,sc,ei,ej,ii,jj,ds_in,type_win);

%this computes only g1111 * s_11
% so for other i and j, it has to be adjusted and all should be added later


%%% write the total_interp out and recon_interp_out to file
% but downsample them first because they were done for the 1024 case.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ds = 1024/128;
windowed_conv = total_interp_out(1:ds:end,1:ds:end,:); 
windowed_model_conv = recon_interp_out(1:ds:end,1:ds:end,:); 


t = squeeze(deltaE(ei,ej,:,:,:)) ;
t = t +windowed_conv; 
deltaE(ei,ej,:,:,:) = t; 

t2 = squeeze(deltaE2(ei,ej,:,:,:)) ;
t2 = t2 +windowed_model_conv;   
deltaE2(ei,ej,:,:,:) = t2;   

end
end


end
 end

toc
%%% currently windowed_conv is only 3d, so repeating data %%
tic
%% for the fortran %%
%for ei = 1:3
%  for ej=1:3
disp('writing')
for kz = 1:nn(3)


   for i=1:3
      for j= 1:3
         i, j
        fname = sprintf('WC%d%d%d.out', (kz + 100) ,  i, j);
        comp = squeeze(deltaE(i,j,:,:,kz));
   

	fileid = fopen(fname, 'w'); 
	[wn wn ] = size(comp);
	for kyy =1:wn
 	 for kxx = 1:wn

     	fprintf(fileid, '%f\n', comp(kxx,kyy));
  	end
	end

	fclose(fileid);
   end
end


 end
%end
%end

toc
quit
