clear all
close all



lnoise1=0.9;	%bpass
lnoise2=0.9;
lobject1=5;    %bpass
lobject2=5;
th1=7.0;	%pkfnd
th2=7.0;
sz1=5;	  %cntrd
sz2=5;
maxdisp=100;	%track
dt=1;   %microseconds
%pixmm=7.11;

x0=80;    
y0=50;	  %pixels
x=1900;
y=1900;

theta=-88.5; 
rect=[x0 y0 x y];




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





k=0;    %loop counter. When k=n, don't clear imageA imageB imgA imgB fobjectA fobjectB (line 227) 
n=length(selectedFiles);


for i=1:2:n 
    
    if i==1        
        image_number=1        
    else        
        image_number=i-image_number         %for more than 1 frames: image_number=i-image_number
    end
    
    flag=1;
    
    imageA(:,:) = imread(selectedFiles{i});
    imageB(:,:) = imread(selectedFiles{i+1});
    
    imageA=imrotate(imageA,theta);
    imageB=imrotate(imageB,theta);



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
    
    
%     %%%%vertical damaged pixel line%%%%    
%     for ipC=1:281
%         
%         for jpC=1:4
%             
%             imageA(11+ipC,1252+jpC)=0;  % (y,x)
%             imageB(11+ipC,1252+jpC)=0;
%                         
%         end
%         
%     end

    
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
    imgB(:,:)=bpass(cimageB(:,:),lnoise2,lobject2);

            
    fobjectA(:,:)=pkfnd(imgA(:,:),th1,sz1);    % finds local maxima in an image to pixel level accuracy. See pkfnd.m for more info.          
    fobjectB(:,:)=pkfnd(imgB(:,:),th2,sz2);
        
%     if i==1
%         figure,
%         imagesc(imgA(:,:));
%         colormap(gray)
%         hold on
%         plot(fobjectA(:,1),fobjectA(:,2),'or')
%         hold off
%         impixelinfo;
%     end
               
    sfobjectA=size(fobjectA(:,:));    %e.g. sfobjectA = 127	 2
    sfobjectA=sfobjectA(1);    %Number of local maxima found with pkfnd e.g. sfobjectA = 127
                
    sfobjectB=size(fobjectB(:,:));
    sfobjectB=sfobjectB(1);
    
    
    if i==1
        size_imgA=size(imgA,1);

        for i=1:size_imgA
            for j=1:size_imgA
              peaks_A(i,j)=0;   %Create a black image same size as imgA   
              peaks_B(i,j)=0;
            end
        end
    end
    

    for i=1:sfobjectA
        
        d=fobjectA(i,1);    %fobjectA is a matrix with 2 columns and sfobjectA rows. It contains the coordinates of local maxima.
        e=fobjectA(i,2);
        peaks_A(e,d) = peaks_A(e,d)+1;    %Thus, we say here e.g. go to peaks_A(1175,1146) and increase intensity by 1.

    end


    for i=1:sfobjectB 
        
       d=fobjectB(i,1);
       e=fobjectB(i,2);
       peaks_B(e,d) = peaks_B(e,d)+1;

    end
  
    
    k=k+1;
    
    if k~=n       
        clear imageA imageB imgA imgB fobjectA fobjectB
    end
          

end



oneA=zeros(size_imgA); %oneA=peaks_A;    %pixels encounteres only once, hopefully only particles
oneB=zeros(size_imgA);
noiseA=zeros(size_imgA);    %the noise that will consist the backround image in main code.
noiseB=zeros(size_imgA);



for i=1:size_imgA
    for j=1:size_imgA
        if peaks_A(i,j)== 1
            oneA(i,j)=500;    %From the peaks_A keep only the local maxima found once (particles)               
        end
        if peaks_B(i,j)== 1
            oneB(i,j)=500;                  
        end
        if peaks_A(i,j)> 1
            noiseA(i,j)=1000;   %From the peaks_A keep only the local maxima founde more than once (noise)
        end        
        if peaks_B(i,j)> 1
            noiseB(i,j)=1000;
        end
    end
end




figure,%('Color',[1 1 1]),
imagesc(oneA);    % Scale data and display as image. imagesc(...) is the same as IMAGE(...) except the data is scaled to use the full colormap.
colormap(gray);    %Color look-up table.
xlabel({'pixels'});    % adds text "pixels" beside the X-axis on the current axis.
ylabel({'pixels'});    % adds text "pixels" beside the Y-axis on the current axis.
set(gca,'DataAspectRatio',[1 1 1]);    %gca: Get handle to current axis. The first element specifying the intensity of red light, the second green, and the third blue.  Color intensity can be specified on the interval 0.0 to 1.0. For example, [0 0 0] is black, [1 1 1] is white
impixelinfo; 


    
figure,%('Color',[1 1 1]),
imagesc(noiseA);    % Scale data and display as image. imagesc(...) is the same as IMAGE(...) except the data is scaled to use the full colormap.
colormap(gray);    %Color look-up table.
xlabel({'pixels'});    % adds text "pixels" beside the X-axis on the current axis.
ylabel({'pixels'});    % adds text "pixels" beside the Y-axis on the current axis.
set(gca,'DataAspectRatio',[1 1 1]);    %gca: Get handle to current axis. The first element specifying the intensity of red light, the second green, and the third blue.  Color intensity can be specified on the interval 0.0 to 1.0. For example, [0 0 0] is black, [1 1 1] is white
impixelinfo;

figure,%('Color',[1 1 1]),
imagesc(noiseB);    % Scale data and display as image. imagesc(...) is the same as IMAGE(...) except the data is scaled to use the full colormap.
colormap(gray);    %Color look-up table.
xlabel({'pixels'});    % adds text "pixels" beside the X-axis on the current axis.
ylabel({'pixels'});    % adds text "pixels" beside the Y-axis on the current axis.
set(gca,'DataAspectRatio',[1 1 1]);    %gca: Get handle to current axis. The first element specifying the intensity of red light, the second green, and the third blue.  Color intensity can be specified on the interval 0.0 to 1.0. For example, [0 0 0] is black, [1 1 1] is white
impixelinfo;    

imwrite(noiseA,'00backroundimageA.tiff');
imwrite(noiseB,'00backroundimageB.tiff');
imwrite(oneA,'particles.tiff');
imwrite(oneB,'particles.tiff');

