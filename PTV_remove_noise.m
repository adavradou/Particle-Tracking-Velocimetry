clear all
close all

% backroundimage_A = imread('00backroundimageA.TIFF');
% backroundimage_B = imread('00backroundimageB.TIFF');

backroundimage_A = imread('00with_bpass_th7_4900_A.TIFF');
backroundimage_B = imread('00with_bpass_th7_4900_B.TIFF');


% backroundimage_A = imread('00with_bpass_th7_A.TIF');
% backroundimage_B = imread('00with_bpass_th7_B.TIF');

% backroundimage_A = imread('00with_bpass_th5_A.TIF');
% backroundimage_B = imread('00with_bpass_th5_B.TIF');

%%%%%%%%%%%%%%%%%%fopen a single file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = uigetdir;
Files = dir(folder);
files = {Files.name};isdir = [Files.isdir];
files(isdir) = [];
s = listdlg('ListString', files);
selectedFiles = files(s);

if iscell(selectedFiles) == 0 %if select only one file, then the data will not be a cell  
    
    image = imread(selectedFiles);
    if image==-1                % If this returns a -1, we did not open the file successfully.
       error('File not found or permission denied');
    end        
        
end
%%%%%%%%%%%%%%%%%%fopen a single file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



lnoise1=0.9;	%bpass
lnoise2=0.9;
lobject1=5;    %bpass
lobject2=5;
th1=10.0;	%pkfnd
th2=10.0;
sz1=5;	  %cntrd
sz2=5;
maxdisp=100;	%track
dt=1;   %microseconds
% pixmm=7.11;

x0=80;    
y0=50;	  %pixels
x=1900;
y=1900;

theta=-88.5; 
rect=[x0 y0 x y];
size_img=size(backroundimage_A,1);



k=0;    %loop counter. When k=n, don't clear imageA imageB imgA imgB fobjectA fobjectB (line 227) 
flag=1;
n=length(selectedFiles);


for i=1:2:n 
    
%     if i==1        
%         image_number=1        
%     else        
%         image_number=i-image_number    %for more than 1 frames: image_number=i-image_number     
%     end
    
    flag=1;
    
    imageA(:,:) = imread(selectedFiles{i});
    imageB(:,:) = imread(selectedFiles{i+1});
    
    imageA=imrotate(imageA,theta);
    imageB=imrotate(imageB,theta);
    
%     if i==1        
%         figure,%('Color',[1 1 1]),
%         imagesc(backroundimage_A);    % Scale data and display as image. imagesc(...) is the same as IMAGE(...) except the data is scaled to use the full colormap.
%         colormap(gray);    %Color look-up table.
%         xlabel({'pixels'});    % adds text "pixels" beside the X-axis on the current axis.
%         ylabel({'pixels'});    % adds text "pixels" beside the Y-axis on the current axis.
%         set(gca,'DataAspectRatio',[1 1 1]);    %gca: Get handle to current axis. The first element specifying the intensity of red light, the second green, and the third blue.  Color intensity can be specified on the interval 0.0 to 1.0. For example, [0 0 0] is black, [1 1 1] is white.
%         impixelinfo;     %creates a pixel information tool in the current figure.  The pixel information tool displays information about the pixel in an image that the cursor is positioned over.
%     end



%%%%%%%%%%%%%%%%%%%remove "bad" pixels%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%

%remember: i refers to y axis and j to x axis.


    %%%%low damaged pixel line%%%%    
    for ipB=1:250
        
        for jpB=1:2048
            
            imageA(1671+ipB,0+jpB)=0;
            imageB(1671+ipB,0+jpB)=0;
            
        end
        
    end    
    
    %%%vertical damaged pixel line%%%%    
    for ipC=1:282
        
        for jpC=1:8
            
            imageA(1171+ipC,1249+jpC)=0;  % (y,x)
            imageB(1171+ipC,1249+jpC)=0;
                        
        end
        
    end    
    
      %%%%damaged spot%%%%      
    for ipC=1:19
        
        for jpC=1:12
            
            imageA(1342+ipC,1562+jpC)=0;  % (y,x)
            imageB(1342+ipC,1562+jpC)=0;
                        
        end
        
    end

    
    %%%%aerolens upper boundary%%%%  
    for ipE=1:170
        
        for jpE=1:740
            
            imageA(610+ipE,0+jpE)=0;
            imageB(610+ipE,0+jpE)=0;
            
            
        end
        
    end


    %%%%aerolens lower boundary%%%%  
    for ipF=1:160
        
        for jpF=1:735
            
            imageA(1400+ipF,0+jpF)=0;
            imageB(1400+ipF,0+jpF)=0;
                        
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%remove "bad" pixels%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


    cimageA(:,:) = imcrop(imageA(:,:),rect);    %crops the image rectangular
    cimageB(:,:) = imcrop(imageB(:,:),rect);
    
  
    imgA(:,:)=bpass(cimageA(:,:),lnoise1,lobject1);    % Implements a real-space bandpass filter. See bpass.m for more info.   
    imgB(:,:)=bpass(cimageB(:,:),lnoise2,lobject2);    % Implements a real-space bandpass filter. See bpass.m for more info.
    
                
    for ii=1:size_img
        for jj=1:size_img
            if backroundimage_A(ii,jj) > 1
                imgA(ii,jj)=0;    %remove the noise(backround image) from imgA      

                imgA(ii,jj-1)=0;    %remove also +/-1 pixels 
                imgA(ii,jj+1)=0;       
                imgA(ii-1,jj)=0;         
                imgA(ii-1,jj-1)=0;
                imgA(ii-1,jj+1)=0; 
                imgA(ii+1,jj)=0;     
                imgA(ii+1,jj-1)=0;
                imgA(ii+1,jj+1)=0;         
                
%                 imgA(ii,jj-2)=0;    %remove also +/-2 pixels 
%                 imgA(ii,jj+2)=0;                          
%                 imgA(ii-1,jj-2)=0;
%                 imgA(ii-1,jj+2)=0;       
%                 imgA(ii+1,jj-2)=0;
%                 imgA(ii+1,jj+2)=0;                   
%                 imgA(ii-2,jj)=0;
%                 imgA(ii-2,jj-1)=0;
% %                 imgA(ii-2,jj-2)=0;
%                 imgA(ii-2,jj+1)=0;
% %                 imgA(ii-2,jj+2)=0;                
%                 imgA(ii+2,jj)=0;
%                 imgA(ii+2,jj-1)=0;
% %                 imgA(ii+2,jj-2)=0;
%                 imgA(ii+2,jj+1)=0;
% %                 imgA(ii+2,jj+2)=0;               
            end
            
            if imgA(ii,jj)>1000    %if something is very bright, remove it
                imgA(ii,jj)=0;    
                imgA(ii,jj-1)=0;
                imgA(ii,jj+1)=0;       
                imgA(ii-1,jj)=0;          
                imgA(ii-1,jj-1)=0;
                imgA(ii-1,jj+1)=0; 
                imgA(ii+1,jj)=0;        
                imgA(ii+1,jj-1)=0;
                imgA(ii+1,jj+1)=0; 
            end
                
            
            
            if backroundimage_B(ii,jj) > 1
                imgB(ii,jj)=0;    %remove the noise(backround image) from imgB   
                
                imgB(ii,jj-1)=0;
                imgB(ii,jj+1)=0;      
                imgB(ii-1,jj)=0;   
                imgB(ii-1,jj-1)=0;
                imgB(ii-1,jj+1)=0; 
                imgB(ii+1,jj)=0;       
                imgB(ii+1,jj-1)=0;
                imgB(ii+1,jj+1)=0;    
                
%                 imgB(ii,jj-2)=0;
%                 imgB(ii,jj+2)=0;                          
%                 imgB(ii-1,jj-2)=0;
%                 imgB(ii-1,jj+2)=0;       
%                 imgB(ii+1,jj-2)=0;
%                 imgB(ii+1,jj+2)=0;                   
%                 imgB(ii-2,jj)=0;
%                 imgB(ii-2,jj-1)=0;
% %                 imgB(ii-2,jj-2)=0;
%                 imgB(ii-2,jj+1)=0;
% %                 imgB(ii-2,jj+2)=0;                
%                 imgB(ii+2,jj)=0;
%                 imgB(ii+2,jj-1)=0;
% %                 imgB(ii+2,jj-2)=0;
%                 imgB(ii+2,jj+1)=0;
% %                 imgB(ii+2,jj+2)=0;
            end
            if imgB(ii,jj)>1000    %if something is very bright, remove it
                imgB(ii,jj)=0;    
                imgB(ii,jj-1)=0;
                imgB(ii,jj+1)=0;       
                imgB(ii-1,jj)=0;          
                imgB(ii-1,jj-1)=0;
                imgB(ii-1,jj+1)=0; 
                imgB(ii+1,jj)=0;        
                imgB(ii+1,jj-1)=0;
                imgB(ii+1,jj+1)=0; 
            end
        end
    end
    
    
      
    fobject=pkfnd(imgA(:,:),th1,sz1);  
    sfobject=size(fobject(:,:));
    sfobject=sfobject(1);
    
    for sp=1:sfobject    %in order to be a particle, it has to have at least one layer of nonzero pixels around its local maxima
        if imgA(fobject(sp,2)-1,fobject(sp,1)) > 0 ...
            & imgA(fobject(sp,2)+1,fobject(sp,1)) > 0 ...
            & imgA(fobject(sp,2),fobject(sp,1)-1) > 0 ... 
            & imgA(fobject(sp,2)-1,fobject(sp,1)-1) > 0 ... 
            & imgA(fobject(sp,2)+1,fobject(sp,1)-1) > 0 ...
            & imgA(fobject(sp,2),fobject(sp,1)+1) > 0 ...
            & imgA(fobject(sp,2)-1,fobject(sp,1)+1) > 0 ...
            & imgA(fobject(sp,2)+1,fobject(sp,1)+1) > 0 ...    %till here is a radius of 1 non-zero particles
            & imgA(fobject(sp,2),fobject(sp,1)-2) > 0 ...
            & imgA(fobject(sp,2)-1,fobject(sp,1)-2) > 0 ...
            & imgA(fobject(sp,2)+1,fobject(sp,1)-2) > 0 ...
            & imgA(fobject(sp,2),fobject(sp,1)+2) > 0 ...
            & imgA(fobject(sp,2)-1,fobject(sp,1)+2) > 0 ...
            & imgA(fobject(sp,2)+1,fobject(sp,1)+2) > 0;    %till here is a radius of 2 non-zero particles  
            
            
        fobject(sp,:)=fobject(sp,:);

                         
        else 
            fobject(sp,:)=0;    
                        
        end   
    end
    
    nonzero=find(fobject(:,1)~=0);
    fobjectA=fobject(nonzero,:);

    
    
    if fobjectA~=0;
        dlmwrite('ImID_20160408_02_aerolens_P20mbar_1pos.txt',selectedFiles{i},'-append','delimiter','\t');    %write image name with particles in frameA, even if untrackable with B        
        sfobjectA=size(fobjectA(:,:));
        sfobjectA=sfobjectA(1);
                
        fcentroid(:,:)=cntrd(imgA(:,:),fobjectA(:,:),sz1);    %calculates the centroid of bright spots to sub-pixel accuracy.                    
        xyzsA(:,:)=horzcat(fcentroid(:,1:2),1*ones([sfobjectA 1]));
                            
    else                    
        flag=0;   
                 
    end
    
    clear fobject sfobject  
    
    
    
    fobject=pkfnd(imgB(:,:),th1,sz1); 
    sfobject=size(fobject(:,:));
    sfobject=sfobject(1);
    
    for sp=1:sfobject    %in order to be a particle, it has to have at least layer of nonzero pixels around it 
        if imgB(fobject(sp,2)-1,fobject(sp,1)) > 0 ...
            & imgB(fobject(sp,2)+1,fobject(sp,1)) > 0 ...
            & imgB(fobject(sp,2),fobject(sp,1)-1) > 0 ... 
            & imgB(fobject(sp,2)-1,fobject(sp,1)-1) > 0 ... 
            & imgB(fobject(sp,2)+1,fobject(sp,1)-1) > 0 ...
            & imgB(fobject(sp,2),fobject(sp,1)+1) > 0 ...
            & imgB(fobject(sp,2)-1,fobject(sp,1)+1) > 0 ...
            & imgB(fobject(sp,2)+1,fobject(sp,1)+1) > 0 ...    %till here is a radius of 1 non-zero particles
            & imgB(fobject(sp,2),fobject(sp,1)-2) > 0 ...
            & imgB(fobject(sp,2)-1,fobject(sp,1)-2) > 0 ...
            & imgB(fobject(sp,2)+1,fobject(sp,1)-2) > 0 ...
            & imgB(fobject(sp,2),fobject(sp,1)+2) > 0 ...
            & imgB(fobject(sp,2)-1,fobject(sp,1)+2) > 0 ...
            & imgB(fobject(sp,2)+1,fobject(sp,1)+2) > 0;    %till here is a radius of 2 non-zero particles         
                  
            
        fobject(sp,:)=fobject(sp,:);

                         
        else 
            fobject(sp,:)=0;
                        
        end   
    end
    
    nonzero=find(fobject(:,1)~=0);
    fobjectB=fobject(nonzero,:);
    sfobjectB=size(fobjectB(:,:));
    sfobjectB=sfobjectB(1);
    
    
    if fobjectB(:,:)~=0;               
                
        sfobjectB=size(fobjectB(:,:));
        sfobjectB=sfobjectB(1);
                   
        fcentroid2(:,:)=cntrd(imgB(:,:),fobjectB(:,:),sz2);    %calculates the centroid of bright spots to sub-pixel accuracy.
        xyzsB(:,:)=horzcat(fcentroid2(:,1:2),2*ones([sfobjectB 1]));                                        
               
    else                
        flag=0;
                %                                 delete(selectedFiles{i},selectedFiles{i+1});                
    end           
    
    
        
    
%     figure,
%     imagesc(imgA(:,:));
%     colormap(gray)
%     hold on
%     plot(fobjectA(:,1),fobjectA(:,2),'or')
%     hold off
%     impixelinfo;
% 
%     
%     
%     figure,
%     imagesc(imgB(:,:));
%     colormap(gray)
%     hold on
%     plot(fobjectB(:,1),fobjectB(:,2),'or')
%     hold off
%     impixelinfo;
                    
    clear fobject sfobject fobjectA sobjectA fobjectB sobjectB  fcentroid fcentroid2    
    
    
    

    

                
           
   
            
    
if flag~=0;    %if the pkfnd has found coordinates of local maxima                 
            
            
            xyzs2(:,:)=vertcat(xyzsA(:,:),xyzsB(:,:));    % Vertical concatenation of matrices xyzsA and xyzsB.
            
            trp(:,:)=track(xyzs2(:,:),maxdisp);    %Constructs n-dimensional trajectories from a scrambled list of particle coordinates determined at discrete times
            
            clear xyzsA xyzsB xyzs2

            
            
            sztrp=size(trp(:,:));       
            sztrp=sztrp(1);        

            if sztrp>max(trp(:,4));                 
                
                for m=2:sztrp
                    
                    if trp(m,4)==trp(m-1,4);         

                        xp(:,m)=(trp(m,1)+trp(m-1,1))/2;
                        yp(:,m)=(trp(m,2)+trp(m-1,2))/2;
                        vx(:,m)=(trp(m,1)-trp(m-1,1))/dt; %*pixmm/dt;
                        vy(:,m)=(trp(m,2)-trp(m-1,2))/dt; %*pixmm/dt;
                        magn(:,m)=sqrt(vx(:,m).^2+vy(:,m).^2);                            
                        
                    else
                        
                        xp(:,m)=0;
                        
                    end
                    
                end
                
                xp(:);
                sxp=size(xp(:));
                sxp=sxp(1);
                
                loc=find(xp(:)~=0);
                
                if loc~=isempty(loc);%                     
                    
                    vecdat(:,1) = xp(:,loc); %xp(:,loc)*(pixmm/1000);
                    vecdat(:,2) = yp(:,loc); %(y-yp(:,loc))*(pixmm/1000);
                    vecdat(:,3) = vx(:,loc);
                    vecdat(:,4) = vy(:,loc);
                    vecdat(:,5) = magn(:,loc);
                    
                
%                 vecscale=0.5;
%                 figure('Color',[1 1 1]),
%                 quiver(vecdat(:,1),vecdat(:,2),vecdat(:,3),vecdat(:,4),vecscale);
%                 xlabel({'x (mm)'});
%                 ylabel({'y (mm)'});
%                 axis([0 1968 0 1998]);
                    
                    k=k+1;
                    
                    if k==1
                        
                        addvec=vecdat;
                        
                    else
                        
                        addvec=vertcat(addvec,vecdat);%                         
                        
                    end
                    
                    
                else
                    
                    flag=2;
                    disp('Warning: Particle Tracking was not possible');
                    %                                 dlmwrite('20130429_case01_no_tracking_possible.txt',selectedFiles{i},'-append','delimiter','\t');
%                       delete(selectedFiles{i},selectedFiles{i+1});
                    
                end
                
            else
                
                flag=2;
                disp('Warning: Particle Tracking was not possible');
                %                         dlmwrite('20130429_case01_no_tracking_possible.txt',selectedFiles{i},'-append','delimiter','\t');
%                   delete(selectedFiles{i},selectedFiles{i+1});
                %             
            end
            
            if flag~=2;
                dlmwrite('201600523_02_aerolens_P20mbar_1pos.txt',addvec,'\t');
            end
        
else
            clear xyzsA xyzsB
            disp('Warning: Particles were not found');
        
%           delete(selectedFiles{i},selectedFiles{i+1});
        
        
    
end
    
    clear trp sztrp vx vy magn xp yp sxp loc vecdat
    
    clear imageA imageB imgA imgB
    
end


