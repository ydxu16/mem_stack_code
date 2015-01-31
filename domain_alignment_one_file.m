%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration plot     %%
%% Nov  2014              %%
%% Yuanda Lawrence Xu     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

set(0,'DefaultAxesFontSize', 14);
set(0,'DefaultTextFontSize', 14);

place_x = 40;
place_y = 80;
ave = 0.4;
var = 0.05;
Nx=64;    %size
Ny=64; 
dy = 0.5;
num_layers = 5;
% Numpics indicate how many output you have, so how many picture 
% you have
Numpics = 80;
Numanaly = 20;
Nstep = 400000;
dt = 0.001;
Nprint = Nstep/Numpics ;
dir = '/Users/ydxu16/Desktop/data/matlab_data_align/';

psi_value = cell(num_layers,1);
% first load all the files, then separate them as times

for i = 1:1:num_layers
    ii = i - 1;
    filename = strcat(dir,'file_membrane_', int2str(ii),'_psi');
    psi_value{i} = load(filename);
end

 for p = 1:1:Numpics+1 % p=1:1:(Numpics+1)
        pp = p - 1;
        for index = 1:1:num_layers
          index_i = index - 1;
          matrix_index = psi_value{index}( ((pp*Nx*Ny+1):((pp+1)*Nx*Ny) ));
 
          for i = 1:1:Nx
            l=i-1;
            for j=1:1:Ny
                % z(x,y,mem_index) at time p
                z(i,j,index)=matrix_index(l*Ny+j,1);        
            end
          end
        end 
    
    % calculate the average(over different position) correlation length for membranes
    for cor_len = 1:1:4
        M_sum = zeros(Nx,Ny);
        for i = 1:1:num_layers-cor_len
            r1 = reshape(z(:,:,i),1,Nx*Ny);
            r2 = reshape(z(:,:,i+cor_len),1,Nx*Ny);
            M_sum = M_sum+( (z(:,:,i)-ave) .* (z(:,:,i+cor_len)-ave))./std2(z(:,:,i))./std2(z(:,:,i+cor_len));  
        end;
            M_sum = M_sum./(num_layers-cor_len);%./(num_layers-cor_len-1);
            ave_cor(cor_len,p) = sum(M_sum(:))/Nx/Ny; 
            
    end;
    
    %calculate average(over different position and membranes) order parameter ^ 2 of all
    %the membranes
    ave_order_mem = zeros(Nx,Ny);
    for i = 1:1:num_layers
        ave_order_mem = ave_order_mem + z(:,:,i).*z(:,:,i);
    end
    ave_order_mem = ave_order_mem ./ num_layers;
    ave_order(p) = sum(ave_order_mem(:))/Nx/Ny;
    
    
    
 %calculate the area and the length of boundary
 for index = 1:1:num_layers
    bw = z(:,:,index);
    level = graythresh(bw);
    bw = im2bw(bw,level);
    bw = ~bw;
    cc = bwconncomp(bw, 4);
    cc.NumObjects;
    graindata = regionprops(cc, 'basic');
    perimeter = regionprops(cc, 'perimeter');
    total_area = 0;
    total_perimeter = 0;
    for index1 = 1:1:cc.NumObjects
         total_area = total_area + graindata(index1).Area;
         total_perimeter = total_perimeter + perimeter(index1).Perimeter;
    end;
    ave_area = total_area/cc.NumObjects;
    domain_size(index,p) = total_area/ total_perimeter;
 end;
     
 end;
 
 %plot the figure
 p = 1:1:Numpics+1; pp = (p-1) * Nstep/Numpics;
 cor_len = 1:1:4;
 hold on;
 plot(pp,ave_cor(cor_len,p));
 xlabel('time');
 ylabel('correlation');
 name1 = strcat('Correlation increases with time in ',int2str(num_layers),' membrane system');
 title(name1);
 name2=strcat('Correlation_time_',int2str(num_layers),'_mem.fig');
 savefig(name2);
 
 hold off;
 
 p = 1:1:Numpics+1;
 pp = (p-1) * Nstep/Numpics;
 plot(pp,ave_order);
 name1 = strcat('AveOrder square increases with time in ', int2str(num_layers), ' membrane system');
 title(name1);
 name2 = strcat('Ave_Order^2_time_',int2str(num_layers),'_mem.fig');
 savefig(name2);
 
 index = 1:1:num_layers;
 p = 1:1:Numpics+1;
 pp = (p-1) * Nstep/Numpics;
 plot(pp,domain_size(index,p));
 xlabel('Time');
 ylabel('Average Domain Size(Total Area/Total length of interface)');
 name1 = strcat('Average Domain Size increases with time in ',int2str(num_layers),' membrane system');
 title(name1);
 name2 = strcat('Average Domain Size in ',int2str(num_layers),'_mem system.fig');
 savefig(name2);
 
%          
%     %calculate the continus phase alignment
%          contin = zeros(num_layers,1);
%          p = 1;
%          r = 1;
%          
%          for index = 1:1:num_layers-1
%           
%            if( z(place_x - 1, place_y - 1, index) * z(place_x - 1, place_y - 1, index+1) > 0)
%                index = index+1;
%                r = r+1;
%               
%            
%            else 
%                 contin(p)=r;
%                 p = p+1;
%                 r=1;  
%            end;
%            
%          end;
%          contin(p) = num_layers;
%          for i_p = 1:1:p-1
%              contin(p)= contin(p) - contin(i_p);
%          end;
%          
     



     %plot the phase at place_x place_y grid for different membranes
    % index = 1:1:num_layers;
    % a=  reshape( z(place_x -1, place_y-1,:) , 1, num_layers);
     % bar(index,a);
        
    %   xname = strcat('Position at (', int2str(place_x),', ',int2str(place_y), ')');
    %   xlabel(xname);
    %   ylabel('Phase');
       

