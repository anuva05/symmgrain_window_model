%works for single fiber examples
%clear
close all
%% extract data from SG files 
N=128; 
n3 = 4;

disp('getting stress')
%foldername='twograin_2d/data';
%plot(micros)
%plot(micros(:,1))

% phi 1 - 1
% phi 2 - 2
% phi 3 - 3
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
% x - 4
% y - 5
% z - 6
%grain id= 7
% phase id = 8






 
    sq_plane = zeros(N,N); %holder for squaring the components
    
     dotc ='.';
     outext = 'out';
    for i =1:3
        for j=1:3
            
            
          
              
          
                for z =1:n3
                        filenum= 100+z;
                        compnum= 10*i + j;
                        full_num = filenum*100 + compnum; 
                        nums=num2str(full_num);
                       filename= sprintf('SG%d%d.out',filenum, compnum)
                       
 %strange but ok 
                     %  filename= sprintf('%s/SG%s%d.%s',foldername, ,4, outext); %strange but ok  
                     tmpsg = importdata(filename);

                        tmpim=reshape(tmpsg,[N,N]);



                        sgCubeComp(i,j,:,:,z)= tmpim;
                        
                    

            
                 end
        end
    end
       

save p128_2grain.mat
