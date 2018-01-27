function [total_interp_out,recon_interp_out] = compress_update_pipeline2d(c0,s0,nn, s11 , ei, ej,ii,jj, ds_in,type_win)
%% function to act as plug in to the iterative matlab solver version
 % data will be large eg. 1024 x 1024  eventually but currently only
 % 128x128 for testing 
% inputs:
% c0,s0 - stiffness and compliance of the material
% nn - size of the microstructure in the format [sx sy sz] 
% s11 - 3d data (one single component of any tensor on the 3d volume)
% ei,ej - refer to the component of strain being computed e_{ij}
% ii  - indexes component of sg being used
% jj - indexes component of sg, same as ii.
% ds_in  - downsample factor desired for grain interior
% set window type
% type_win=1 ; %tukey win
% type_win=2; % trap win
% NOT YET SET - window parametere inputs
% outputs: 
%total_interp_out = putting together grain interior and exterior that have
%been convolved separately. no modeling
% recon_interp_out = grain interior has been modeled and compressed, then
% reconstructed. then, interpolated and added to surface convo. so bound to
% have more errors

nn = [1024 1024 4] ;
n1= nn(1);
n2=nn(2);
n3=nn(3);
N = nn(1);


%% generate data by interpolating 128 data 
%% REMOVE THIS PART IN FINAL VERSION 
ds_fake = n1/128;
green_sz = 128;%size of greens core
[X Y ] = meshgrid(-n1/2:n1/2-1 , -n2/2:n2/2-1);
[XIN YIN] =meshgrid(1:nn(1), 1:nn(2));
[XDS YDS] =meshgrid(1:ds_fake:nn(1), 1:ds_fake:nn(2));





%% load data. TEMPORARY FAKE DATA CREATION 

for z=  1:n3
  temp = squeeze(s11(:,:,1));
  strain(:,:,z) =interp2(XDS, YDS, squeeze(s11(:,:,z)), XIN, YIN );%,'linear',4); 
end

%% temporary solution to nan
strain(:,1018:1024,:) = strain(:,1011:1017,:);
strain(1018:1024,:,:) = strain(1011:1017,:,:);


%%%%%%%%%%%% get surface and interior using Ewald window %%%%%%%%%%%%%%%%%
N=n1/2; 

wtemp = zeros(N,N);
[x ] = meshgrid(-N/2: (N/2 -1 ));

splane = strain(:,:,1); 

b1=20;
b2 = 40;

if(type_win ==1)
%% using Tukey win
r= 0.2; 
wc= tukeywin(n1/2,r);
wr=tukeywin(n1/2,r);
[maskr,maskc]=meshgrid(wr,wc);
wtemp=maskr.*maskc;

elseif(type_win==2)
wtemp= trapmf(x , [(-N/2+  b1) (-N/2 + b2)  (N/2 - b2)  (N/2 - b1) ]); %offset should be less than radius

wtemp = wtemp.*wtemp';

 end
    
rect_win = zeros(n1/2,n1/2);
rect_win(b1+1:end-b1,b1+1:end-b1) = 1 ;

rect_win_3d = repmat(rect_win,[1 1 n3]);
int_win_2d = repmat(wtemp,[2 2]);
rect_win_all_2d = repmat(rect_win,[2,2]);


ext_win_2d = 1 - int_win_2d; 

%%%%%%%%%%%%%%% get interior and boundary %%%%%%%%%%%%%%%%%


%%3d
for z = 1:nn(3)
    splane = strain(:,:,z);

    interior3d(:,:,z) =splane.*int_win_2d;
    
    interior_rect3d(:,:,z) = splane.*rect_win_all_2d; 
    exterior3d(:,:,z) = splane.*ext_win_2d; 
    
    int_win_3d(:,:,z)=int_win_2d;
   
    
end
    
%% here we simulate the modeling and transmission of interior_rect3d

%%% extract each grain interior in the rect window area and 
 g1out=( 1 -rect_win_3d).*strain(1:n1/2,1:n1/2,:);
 g2out= ( 1 -rect_win_3d).*strain(1+n1/2:end,1+n1/2:end,:);
 g3out= ( 1 -rect_win_3d).*strain(1:n1/2,1+n1/2:end,:);
 g4out=( 1 -rect_win_3d).*strain(1+n1/2:end,1:n1/2,:);
g1= strain(1:n1/2,1:n1/2,:);
g2= strain(1+n1/2:end,1+n1/2:end,:);
g3= strain(1:n1/2,1+n1/2:end,:);
g4= strain(1+n1/2:end,1:n1/2,:);

%equivalent to recwin but takes off the edges
g1in = g1(b1+1:end-b1,b1+1:end-b1,:);
g2in = g2(b1+1:end-b1,b1+1:end-b1,:);
g3in = g3(b1+1:end-b1,b1+1:end-b1,:);
g4in = g4(b1+1:end-b1,b1+1:end-b1,:);

%% modeling efforts
 [grecon, params] = model_before_tx(g1in); %generates model for 2D plane
 % can do similar for other grains
% these model parameters are transmitted
%% generate greens func



matn = sprintf('Gfft%d%d%d%d.mat',ei,ej,ii,jj);

Gmat=  load(matn);
Gfft = Gmat.G;


%% PERFORM convolutions. single component 

%convo_ds = zeros(n1/ds ,n2/ds, n3 );
i =ei;
j  = ej ;
k= ii;
l=jj;

    greens = ifftshift(ifftn((Gfft(:,:,:))));
    greens_fft_1111 = squeeze(Gfft(:,:,:));
    
    greens_ds_in = greens(1:ds_in:end, 1:ds_in:end, :);     
       
  %  greens_core = greens(18:49, 18:49, :); %greens
     greens_core= greens(nn(1)/2 - (green_sz/2- 1):nn(1)/2 + green_sz/2,nn(2)/2 - (green_sz/2- 1):nn(2)/2 + green_sz/2, : ); %greens 
    %18:49


%% full res with surface
convo_surf =( ifftn(fftn(exterior3d).*fftshift(fftn(greens)))); %this does not
%give border results as expectd
fs=fftshift(convo_surf);

%%%%%%%%%%%%%%%%%%%%% downsampled conv for interior  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
grain_field_ds = interior3d(1:ds_in:end, 1:ds_in:end,:);


%% grain
convo_ds_alt =( ifftn(fftn(grain_field_ds).*fftshift(fftn(greens_ds_in)))); %this does not
%give border results as expectd
f= fftshift(convo_ds_alt);



%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% 
%% actual  

scomp_fft = fftn(strain).*fftshift(greens_fft_1111); %field
actual_conv = ifftn(scomp_fft);

%% interp grain 


%% do  fourier interp 
f1 = interpft(f,n1);
f2 = interpft(f1,n2,2);
f3 = interpft(f2,n3,3);


%% choose which type of interpolation. fourier interp gives better results
%total_interp_out = interpc +  fs; 
total_interp_out = f3 +  fs;


%% for putting back together simulation, we should re-model the interpolated conv. Or model it in downsampled form, transmit, then recreate. 
 % do this for each grain

 fp1 = f(:,:,1); %SHOULD BE FOR EVERY GRAIN.
 %model the central part
  csz= size(fp1,1);
 win_ds = wtemp(1:ds_in:end,1:ds_in:end);
 unwin = win_ds; 
 unwin(win_ds == 0) =1 ; 
 

 %% PICKING THE GRAIN
 recon_grain = zeros(size(fp1)); 
 for g = 1:4
    gr = mod(g,2);
    qx= (1-gr)*(csz/2)+1:csz/2 + (1-gr)*(csz/2) ;
    qy = (1-gr)*(csz/2)+1:csz/2 + (1-gr)*(csz/2) ;
    
    
      if (g/2 <= 1 )
    qx= 1:(csz/2) ;
    qy = (1-gr)*(csz/2)+1:csz/2 + (1-gr)*(csz/2) ;
    
    end
    
    if(g/2 > 1)
     qx=(csz/2)+1 : csz; ;
    qy = (1-gr)*(csz/2)+1:csz/2 + (1-gr)*(csz/2) ; 
        
    end
    
    
   
    gr1c= fp1(qx,qy);
 
    %%% unwindow  %%%

     %gr1c= f(1:csz/2,1:csz/2,1:csz3/2 ) ;  %grain1 after convolution 
     unwindowed_grain = gr1c./unwin;  % "unwindow" ie restore the edges for ease of modeling
     
    % temp=( rect_win(1:ds_in:end,1:ds_in:end).*unwindowed_grain);
    % r=temp(abs(temp)>10^-23);
     rectds = rect_win(1:ds_in:end,1:ds_in:end); 
     r =  unwindowed_grain(rectds==1); %FIND A WAY TO REMOVE THE OUTSIDE EDGES
     %since sq grains
      lr = length(r);
      res = sqrt(lr);
     
      %grain_central = reshape(r,res,res);
     grain_central = unwindowed_grain(7:123 , 7:123); %only for ds = 4
      [recon_post, params_post] = model_after_conv(grain_central); 
        diff = recon_post - grain_central; 
       norm_err_post_conv(g) = norm(diff(:))/norm(grain_central(:))
       put_back = unwindowed_grain;
       put_back(7:123,7:123) = recon_post;
       recon_grain(qx,qy) = put_back.*unwin ; %re-multiply with the window edges
 end


recon_f = repmat(recon_grain,[ 1 1 n3]); %because 2d, we only need to model one plane. all other planes have the same values of stress/strain
% do  fourier interp 
rf1 = interpft(recon_f,n1);
rf2 = interpft(rf1,n2,2);
rf3 = interpft(rf2,n3,3);


%% choose which type of interpolation. fourier interp gives better results
%total_interp_out = interpc +  fs; 
recon_interp_out = rf3 +  fs;
