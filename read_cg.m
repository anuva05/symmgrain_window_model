nn = [128 128 4];
fid = fopen('cg.out','r'); %# open new csv file
% get only single cg line from the file corresponding to tempvec
% 1, 18, 35 etc
cgpart = zeros(1,81);

avgcg = zeros(1,81); 


while  ~feof(fid)
%% now read line for grain %%
    line = fgets(fid);
    comp = line(2:3);
    if strcmp(comp,'cg')

    temp=strsplit(line,'=');
%

   ctr = 1;
    while ctr<17
%     %second line is data=
   line = fgets(fid);
     r = strsplit(line, ' ');
     vals = str2double(r);
     cgpart(5*(ctr -1 ) + 1: 5*(ctr)) = vals(2:6);
     ctr = ctr + 1;
    end
  %17th line is tricky 
   line = fgets(fid);
    r = strsplit(line, ' ');
    vals = str2double(r);
   cgpart(end) = vals(2);
end

avgcg = avgcg + cgpart ;

end


avgcg = avgcg/(nn(1)*nn(2)*nn(3)) ;
fclose(fid);

save avgcg.mat avgcg
