function [recon, params] = model_after_conv(gc)

gp1 = gc;


grsz = size(gp1,1) ;

xcof = 1:size(gp1,1);
ycof= 1:size(gp1,2);

xcof  =xcof - grsz/2; 
ycof = ycof - grsz/2;

theta = 45; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% Rotate your point(s)

cor = [xcof ; ycof]; 

for i = 1:size(cor,1)
    
    corrot(:,i) = R*cor(:,i) ;
end


%[X Y] =meshgrid(corrot(1,:),corrot(2,:));
%
[X Y] =meshgrid(xcof,ycof);


for t = 1:size(gp1,1)
    
    curve = gp1(:,t); 
    xcof = 1:size(gp1,1); 
    xcof = xcof - grsz/2; 
    [mod ]= fit(xcof', curve, 'poly4'); 
    coef_val(t,:) = coeffvalues(mod);
end
    

% get models for the parameters


for u=1:5
    p = coef_val(:,u); 
    curve = p; 
    if u==4
      mod3 = fit(xcof', curve, 'poly5'); 
    end
    
    if u ==1 
        mod0= fit(xcof', curve, 'poly4'); 
    end
    
      if u ==2
        mod1 = fit(xcof', curve, 'poly4'); 
      end
    
      if u ==3
        mod2 = fit(xcof', curve, 'poly4'); 
      end
    
      if u == 5
        mod4 = fit(xcof', curve, 'poly4'); 
    end
    
    
    
    %p0,p1,p2 = poly 4
    % p3 = poly5 
    % p4 = poly4 
end

% 
    xcof = 1:size(gp1,1); 
 xcof = (xcof - grsz/2); 
     X =xcof; 
     ycof = 1:size(gp1,2); 
 ycof = ycof - grsz/2 ; 

recon=zeros(size(gp1)); 
for yy = 1:length(ycof)
    
    yvar = ycof(yy); 
    
    %get value of params at this y
    p0 = feval(mod0,yvar);
    p1 = feval(mod1,yvar);
     p2= feval(mod2,yvar);
   p3 = feval(mod3,yvar);
   p4 = feval(mod4,yvar);
                
          % evaluate x for this set of params      
 recon(:,yy) = p0.*X.^4 + p1.*X.^3 + p2.*X.^2 + p3.*X + p4;
                
end

params{1} = mod; % family of curves eq
%models for parameters
params{2} =mod0;
params{3} =mod1; 
params{4} =mod2; 
params{5} =mod3; 
params{6} =mod4; 