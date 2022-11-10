clc
close all
clear all
warning off
%% inisalizing variables
%images need to be in one folder
% full directory to folder
dirs_b={'G:\Shared drives\Rao Lab\Stem Cell Lab\3D paper\10-30-22\Decidualized stromal rep 1 DAPI_1'};
dirs_g={'G:\Shared drives\Rao Lab\Stem Cell Lab\3D paper\10-30-22\Decidualized stromal rep 1 HLA-G_1'};
dirs_r={'G:\Shared drives\Rao Lab\Stem Cell Lab\3D paper\10-30-22\Decidualized stromal rep 1 Vimentin_1'};
conditions={'1'};
 xlxfilename='IC'; % name/location for genrated excel sheet will apper in the folder matlab is open too
 folder='G:\Shared drives\Rao Lab\Stem Cell Lab\Code';
 if ~exist(folder, 'dir')
     mkdir(folder);
 end
 sheetname=fullfile(folder, xlxfilename);
between_slice=10;%in micrometers
% for bpass method you may need to play with these values depending on image taking procedure
%currently set up for high quality main campus images
noise=1;
object=15;
thresh=.05;
totaldepth=800;
pixel_to_micron=1;
p2mt=2048;
%alpha is need for writing to excel file
alpha={'A','B','C','D', 'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
stats={'Average','STDEV','Cell Count','SEM','MAX','T-test'};
% xlswrite(sheetname,stats(1),'Summary',"A2:A2")
% xlswrite(sheetname,stats(2),'Summary',"A3:A3")
% xlswrite(sheetname,stats(3),'Summary',"A4:A4")
% xlswrite(sheetname,stats(4),'Summary',"A5:A5")
% xlswrite(sheetname,stats(5),'Summary',"A6:A6")
%%
for condition=1:length(conditions)
    %clear all_pos_data length_all_pos image_list
    dirname=string(dirs_b(condition)+"\");
    image_list=dir(fullfile(dirname,'*.tif'));% Identifies files that end with .tif within that folder
    image_list_size=size(image_list);
    length_all_pos=zeros(image_list_size(:,1),1);
    %% goes through each position
full_size_BW=zeros(6144,6144);
full_size_unprocessed=zeros(6144,6144);
full_size_unprocessed_g=zeros(6144,6144);
full_size_unprocessed_r=zeros(6144,6144);
full_size_I=zeros(6144,6144);
    for position=1:length(image_list)
        clearvars -except full_size_BW_inv image_data full_size_BW full_size_I full_size_unprocessed full_size_unprocessed_g full_size_unprocessed_r Green_new_tot Red_new_tot expresh_new_tot cells_new_tot avgdist cells_z_traveledi stats dirs_b dirs_g dirs_r sheetname position dirname totaldepth image_list noise object thresh between_slice all_pos_data length_all_pos alpha conditions condition all_conditions
        %% initializing the z-stack images
   %% initializing the z-stack images
    fprintf('Reading Image Stack %1.0f\n',position)
   % 'Reading Image Stack ',num2str(position))
    file_s=strcat(dirname,image_list(position).name); % creates a variable with the full path to image # 1 (if there was a 2 inside the parenthesis, that would be image #2)
    file=file_s{1};
    info=imfinfo(file);
    num_images=numel(info);
    TifLink=Tiff(file,'r');
    for i=1:num_images
        TifLink.setDirectory(i);
        FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
    img_3d=FinalImage;
  
    z_max=num_images*between_slice;
    %% preparing and processing the max image
    [max_image(:,:,1),Imax(:,:,1)]=max(img_3d,[],3);
    num_z_stack=size(img_3d);
    
    gray_image=mat2gray(max_image,[0 2^16-1]); % normalizes the image from 0 to 1, the limits should be consistent with the bit depth of the image (2^12-1, in this case since 12 bit)
    
    %% processing max image
    distances_size = bpass(gray_image,noise,object,thresh);
    BW2=imdilate(distances_size,strel('sphere',2));
    BW3=bwareaopen(BW2,100);
    BW4=imfill(BW3,'holes');
    BW5=imclearborder(BW4);
    %%
unprocessed_image=max_image;
if position==1
full_size_BW(1:2048,2049:4096)=BW5(1:2048,1:2048);
full_size_unprocessed(1:2048,2049:4096)=unprocessed_image(1:2048,1:2048);
full_size_I(1:2048,2049:4096)=Imax(1:2048,1:2048);
elseif position==2
full_size_BW(2049:4096,1:2048)=BW5(1:2048,1:2048);
full_size_unprocessed(2049:4096,1:2048)=unprocessed_image(1:2048,1:2048);
full_size_I(2049:4096,1:2048)=Imax(1:2048,1:2048);
elseif position==3
full_size_BW(2049:4096,2049:4096)=BW5(1:2048,1:2048);
full_size_unprocessed(2049:4096,2049:4096)=unprocessed_image(1:2048,1:2048);
full_size_I(2049:4096,2049:4096)=Imax(1:2048,1:2048);
elseif position==4
full_size_BW(2049:4096,4097:6144)=BW5(1:2048,1:2048);
full_size_unprocessed(2049:4096,4097:6144)=unprocessed_image(1:2048,1:2048); 
full_size_I(2049:4096,4097:6144)=Imax(1:2048,1:2048);
elseif position==5
full_size_BW(4097:6144,2049:4096)=BW5(1:2048,1:2048);
full_size_unprocessed(4097:6144,2049:4096)=unprocessed_image(1:2048,1:2048); 
full_size_I(4097:6144,2049:4096)=Imax(1:2048,1:2048);
end
CC=bwconncomp(full_size_BW);
%B=imresize(unprocessed_image, 0.25);
%imshow(B)
%C=rescale(B);%not BW image
%figure(1),imshowpair(B,C,'montage')
%im=C(:,:,:) %BW1 is black when used this
%[L A Bl] = imsplit(rgb2lab(C))
%blueChannel = C;%this is a BW picture
%figure(1),imshow(blueChannel)
%figure(2),histogram(blueChannel)
%blueChannel2=blueChannel>0.3; 
%figure(3),imshow(blueChannel2) %why this is BW picture?
%BW1n=mat2gray(blueChannel2); %removed [0 2^8-1]) and it shows an area now %Q-can I show both L and b in Lab color space
%BW2n=imdilate(BW1n,strel('sphere',1)); % dilates objects
%BW3n=imbinarize(BW1n);
%CC=bwconncomp(BW2);
%BW4n = bwareafilt(BW3n,[1 50]);
    %figure(position);imshow(BW5,'InitialMagnification','fit')
    end
    %%
image_data=regionprops(CC,full_size_unprocessed,'all');
    for j=1:numel(image_data)
        ctrs(j,1:2)=image_data(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
    end
    
    %for position=1:length(image_list)
center=[3072 3072];

 
%%
centerx=mean(ctrs(:,1));
centery=mean(ctrs(:,2));
centroid=[centerx centery];
geocenterx=(max(ctrs(:,1))+min(ctrs(:,1)))/2;
geocentery=(max(ctrs(:,2))+min(ctrs(:,2)))/2;
geocenter=[geocenterx geocentery];
    for j=1:numel(image_data)
        d1(j,1)=pdist([ctrs(j,:);geocenter(1,:)]); % calculates distance (euclidean distance) between centroids of Blue and Color objects (red or green)
    end
    ctrsnew=ctrs;
    d2=d1;
for j=1:numel(image_data)
    if d1(j,1)>2000
        ctrsnew(j,1:2)=nan();
        full_size_BW(image_data(j).PixelIdxList)=0;
        d2(j,1)=nan();
    end
end
d2(isnan(d2))=[];
full_size_BW_inv=zeros(6144,6144);
for j=1:6144
    for k=1:6144
        if full_size_BW(j,k)==0
             full_size_BW_inv(j,k)=1;
        else
     full_size_BW_inv(j,k)=0;
        end
    end
end

%     for j=1:numel(image_data)
%         Area(j,1)=image_data(j).Area; % loop goes thorugh each element and extracts centroid from structure B
%     end
% for j=1:numel(image_data)
%     if Area(j,1)<150
%         ctrsnew(j,1:2)=nan();
%         full_size_BW(image_data(j).PixelIdxList)=0;
%     end
% end
%figure(1)=scatter(ctrsnew(:,1),ctrsnew(:,2));
CCnew=bwconncomp(full_size_BW);
CCback=bwconncomp(full_size_BW_inv);
image_data_b=regionprops(CCnew,full_size_unprocessed,'all');
image_data_back=regionprops(CCback,full_size_unprocessed,'all');
%%
for position=1:length(image_list)

dirname_g=string(dirs_g(condition)+"\");
    image_list_g=dir(fullfile(dirname_g,'*.tif'));% Identifies files that end with .tif within that folder
    image_list_size_g=size(image_list_g);
    length_all_pos_g=zeros(image_list_size_g(:,1),1);
    position_g=position;
    %% goes through each position
    %for position_g=1:length(image_list_g)
        %clearvars -except expresh_new_tot cells_new_tot avgdist cells_z_traveledi stats dirs sheetname position dirname totaldepth image_list noise object thresh between_slice all_pos_data length_all_pos alpha conditions condition all_conditions
        %% inslizing the z-stack images
   %% inslizing the z-stack images
    fprintf('Reading Image Stack %1.0f\n',position_g)
   % 'Reading Image Stack ',num2str(position))
    file_s_g=strcat(dirname_g,image_list_g(position_g).name); % creates a variable with the full path to image # 1 (if there was a 2 inside the parenthesis, that would be image #2)
    file_g=file_s_g{1};
    info_g=imfinfo(file_g);
    num_images_g=numel(info_g);
    TifLink=Tiff(file_g,'r');
    for i=1:num_images_g
        TifLink.setDirectory(i);
        FinalImage_g(:,:,i)=TifLink.read();
    end
    TifLink.close();
    img_3d_g=FinalImage_g;
    [max_image_g(:,:,1),Imax_g(:,:,1)]=max(img_3d_g,[],3);
    num_z_stack_g=size(img_3d_g);
    
    gray_image_g=mat2gray(max_image_g,[0 2^16-1]); % normalizes the image from 0 to 1, the limits should be consistent with the bit depth of the image (2^12-1, in this case since 12 bit)
    
%%
unprocessed_image_g=max_image_g;
if position==1
full_size_unprocessed_g(1:2048,2049:4096)=unprocessed_image_g(1:2048,1:2048);
elseif position==2
full_size_unprocessed_g(2049:4096,1:2048)=unprocessed_image_g(1:2048,1:2048);
elseif position==3
full_size_unprocessed_g(2049:4096,2049:4096)=unprocessed_image_g(1:2048,1:2048);
elseif position==4
full_size_unprocessed_g(2049:4096,4097:6144)=unprocessed_image_g(1:2048,1:2048); 
elseif position==5
full_size_unprocessed_g(4097:6144,2049:4096)=unprocessed_image_g(1:2048,1:2048);  
end
%B_g=imresize(unprocessed_image_g, 0.25);
%imshow(B)
%C_g=rescale(B_g);%not BW image
%figure(1),imshowpair(B,C,'montage')
%im=C(:,:,:) %BW1 is black when used this
%[L A Bl] = imsplit(rgb2lab(C))
%ColorChannel = C_g;%this is a BW picture
%figure(1),imshow(blueChannel)
%figure(2),histogram(blueChannel)
%ColorChannel2=ColorChannel>0; 
%figure(3),imshow(blueChannel2) %why this is BW picture?
%BW1n_g=mat2gray(ColorChannel); %removed [0 2^8-1]) and it shows an area now %Q-can I show both L and b in Lab color space
%BW2n=imdilate(BW1n,strel('sphere',1)); % dilates objects
%BW3n_g=imbinarize(BW1n_g);
%CC=bwconncomp(BW2);
%BW4n = bwareafilt(BW3n,[1 50]);
    %figure(position);imshow(BW5,'InitialMagnification','fit')
    %%
end
    image_data_g=regionprops(CCnew,full_size_unprocessed_g,'all');
    image_data_g_back=regionprops(CCback,full_size_unprocessed_g,'all');
 for position=1:length(image_list)
dirname_r=string(dirs_r(condition)+"\");
    image_list_r=dir(fullfile(dirname_r,'*.tif'));% Identifies files that end with .tif within that folder
    image_list_size_r=size(image_list_r);
    length_all_pos_r=zeros(image_list_size_r(:,1),1);
    position_r=position;
    %% goes through each position
    %for position_r=1:length(image_list_r)
        %clearvars -except expresh_new_tot cells_new_tot avgdist cells_z_traveledi stats dirs sheetname position dirname totaldepth image_list noise object thresh between_slice all_pos_data length_all_pos alpha conditions condition all_conditions
        %% inslizing the z-stack images
   %% inslizing the z-stack images
    fprintf('Reading Image Stack %1.0f\n',position_r)
   % 'Reading Image Stack ',num2str(position))
    file_s_r=strcat(dirname_r,image_list_r(position_r).name); % creates a variable with the full path to image # 1 (if there was a 2 inside the parenthesis, that would be image #2)
    file_r=file_s_r{1};
    info_r=imfinfo(file_r);
    num_images_r=numel(info_r);
    TifLink=Tiff(file_r,'r');
    for i=1:num_images_r
        TifLink.setDirectory(i);
        FinalImage_r(:,:,i)=TifLink.read();
    end
    TifLink.close();
    img_3d_r=FinalImage_r;
    [max_image_r(:,:,1),Imax_r(:,:,1)]=max(img_3d_r,[],3);
    num_z_stack_r=size(img_3d_r);
    
    gray_image_r=mat2gray(max_image_r,[0 2^16-1]); % normalizes the image from 0 to 1, the limits should be consistent with the bit depth of the image (2^12-1, in this case since 12 bit)
    
%%
unprocessed_image_r=max_image_r;
if position==1
full_size_unprocessed_r(1:2048,2049:4096)=unprocessed_image_r(1:2048,1:2048);
elseif position==2
full_size_unprocessed_r(2049:4096,1:2048)=unprocessed_image_r(1:2048,1:2048);
elseif position==3
full_size_unprocessed_r(2049:4096,2049:4096)=unprocessed_image_r(1:2048,1:2048);
elseif position==4
full_size_unprocessed_r(2049:4096,4097:6144)=unprocessed_image_r(1:2048,1:2048); 
elseif position==5
full_size_unprocessed_r(4097:6144,2049:4096)=unprocessed_image_r(1:2048,1:2048);  
end
%B_r=imresize(unprocessed_image_r, 0.25);
%imshow(B)
%C_r=rescale(B_r);%not BW image
%figure(1),imshowpair(B,C,'montage')
%im=C(:,:,:) %BW1 is black when used this
%[L A Bl] = imsplit(rgb2lab(C))
%ColorChannel = C_r;%this is a BW picture
%figure(1),imshow(blueChannel)
%figure(2),histogram(blueChannel)
%ColorChannel2=ColorChannel>0; 
%figure(3),imshow(blueChannel2) %why this is BW picture?
%BW1n_r=mat2gray(ColorChannel); %removed [0 2^8-1]) and it shows an area now %Q-can I show both L and b in Lab color space
%BW2n=imdilate(BW1n,strel('sphere',1)); % dilates objects
%BW3n_r=imbinarize(BW1n_r);
%CC=bwconncomp(BW2);
%BW4n = bwareafilt(BW3n,[1 50]);
    %figure(position);imshow(BW5,'InitialMagnification','fit')
    %%
end
    image_data_r=regionprops(CCnew,full_size_unprocessed_r,'all');
    image_data_r_back=regionprops(CCback,full_size_unprocessed_r,'all');   
 if numel(image_data_b)~=0
    %% processing the max images for comparson to z-stack
%     for j=1:numel(image_data)
%         ctrs(j,1:2)=image_data(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
%     end
%     
%     for j=1:numel(image_data)
%         bbx(j,1:4)=image_data(j).BoundingBox;
%     end
%     for j=1:numel(image_data_g)
%         pxmean(j)=image_data_g(j).MeanIntensity; % loop goes thorugh each element and extracts centroid from structure B
%     end
    %% processing the max images for comparson to z-stack
    Blueim_max=full_size_BW;
    Blue_max=image_data_b;
    max_size=size(Blueim_max);
    Imax_z=zeros(max_size(:,1),max_size(:,2));
 %   Px_max=zeros(max_size(:,1),max_size(:,2));
    for j=1:max_size(:,1)
        for k=1:max_size(:,2)
            if Blueim_max(j,k)==1
                Imax_z(j,k)=full_size_I(j,k);
                %Px_max(j,k)=max_image_g(j,k);
            %else Imax_z(j,k)=0;
            end
        end
    end
    
    
%     for j=1:numel(Blue_max)
%         ctrs(j,1:2)=Blue_max(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
%     end
%     
%     for j=1:numel(Blue_max)
%         bbx(j,1:4)=Blue_max(j).BoundingBox;
%     end
%     for j=1:numel(image_data_g)
%         pxmean(j)=image_data_g(j).MeanIntensity; % loop goes thorugh each element and extracts centroid from structure B
%     end
     for j=1:numel(Blue_max)
         pxlist=Blue_max(j).PixelIdxList; % loop goes thorugh each element and extracts centroid from structure B
         k=length(pxlist);
         pxlist_z(1:k,j)=pxlist;
     end
    %%
 %   pxlist_z(pxlist_z==0)=nan();
    Imax_z_filt=reshape(Imax_z,[],1);
    %Px_max_filt=reshape(Px_max,[],1);
    size_max=size(pxlist_z);
    for j=1:size_max(:,1)
        for k=1:size_max(:,2)
            pos=pxlist_z(j,k);
            if pos~=0
            Max_Ind(j,k)=Imax_z_filt(pos);
            %PX_max(j,k)=Px_max_filt(pos);
            else Max_Ind(j,k)=nan();
                %PX_max(j,k)=nan();
            end
        end
    end
    
    Max_avg=mean(Max_Ind,1,'omitnan');
    %MI_avg=mean(PX_max,1,'omitnan');
    cells_new=Max_avg*between_slice;
   % expresh_new=MI_avg;
    cells_z_new=z_max-cells_new;
    
    cells_new_sh=reshape(cells_new,[],1);
   % expresh_new_sh=reshape(expresh_new,[],1);
    s=numel(image_data_b);
    
    for k=1:s
    cells_new_tot(k,position)=cells_new_sh(k,1);
   % expresh_new_tot(k,position)=expresh_new_sh(k,1);
    end
    cells_z_new_tot=z_max-cells_new_tot;
    for j=1:numel(image_data_g)
        Exp_g(j)=image_data_g(j).MeanIntensity; % identifies the particular pixel ID for each object. Identifies where the cells are
    end
    Exp_new_sh=reshape(Exp_g,[],1);
    
    for k=1:numel(image_data_g)
    Green_new_tot(k,position_g)=Exp_new_sh(k,1);
    end
    
    for j=1:numel(image_data_r)
        Exp_r(j)=image_data_r(j).MeanIntensity; % identifies the particular pixel ID for each object. Identifies where the cells are
    end
    Red_new_sh=reshape(Exp_r,[],1);
    
    for k=1:numel(image_data_r)
    Red_new_tot(k,position_r)=Red_new_sh(k,1);
    end
    
    for j=1:numel(image_data_g_back)
        Exp_g_back(j)=image_data_g_back(j).MeanIntensity; % identifies the particular pixel ID for each object. Identifies where the cells are
    end
    Exp_back_sh=reshape(Exp_g_back,[],1);
    
    for k=1:numel(image_data_g_back)
    Green_back_tot(k,position_g)=Exp_back_sh(k,1);
    end
    
    for j=1:numel(image_data_r_back)
        Exp_r_back(j)=image_data_r_back(j).MeanIntensity; % identifies the particular pixel ID for each object. Identifies where the cells are
    end
    Red_back_sh=reshape(Exp_r_back,[],1);
    
    for k=1:numel(image_data_r_back)
    Red_back_tot(k,position_r)=Red_back_sh(k,1);
    end
    %cells_z_new_tot=z_max-cells_new_tot;
else cells_new_tot(1,position)=nan();
     Green_new_tot(1,position_g)=nan();
     Red_new_tot(1,position_r)=nan();
     Green_back_tot(1,position_g)=nan();
     Red_back_tot(1,position_r)=nan();
 end

    Green_new_tot(Green_new_tot==0)=[];
    Red_new_tot(Red_new_tot==0)=[];
    cells_new_tot(cells_new_tot==0)=[];
    Green_back_tot(Green_back_tot==0)=[];
    Red_back_tot(Red_back_tot==0)=[];
    Z_depth_data=reshape(cells_new_tot,[],1);
    Green_data=reshape(Green_new_tot,[],1);
    Red_data=reshape(Red_new_tot,[],1);
    Green_background=reshape(Green_back_tot,[],1);
    Red_background=reshape(Red_back_tot,[],1);
    Red_final=Red_data-Red_background;
    Green_final=Green_data-Green_background;
 end
%%
% 
%  if numel(image_data)~=0
%     %% processing the max images for comparson to z-stack
% %     for j=1:numel(image_data)
% %         ctrs(j,1:2)=image_data(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
% %     end
% %     
% %     for j=1:numel(image_data)
% %         bbx(j,1:4)=image_data(j).BoundingBox;
% %     end
%     for j=1:numel(image_data_r)
%         pxmean(j)=image_data_r(j).MeanIntensity; % loop goes thorugh each element and extracts centroid from structure B
%     end
%     %% processing the max images for comparson to z-stack
%     Blueim_max=BW5;
%     Blue_max=image_data;
%     max_size=size(Blueim_max);
%     Imax_z=zeros(max_size(:,1),max_size(:,2));
%     Px_max=zeros(max_size(:,1),max_size(:,2));
%     for j=1:max_size(:,1)
%         for k=1:max_size(:,2)
%             if Blueim_max(j,k)==1
%                 Imax_z(j,k)=Imax(j,k);
%                 Px_max(j,k)=max_image_r(j,k);
%             %else Imax_z(j,k)=0;
%             end
%         end
%     end
    
    
%     for j=1:numel(Blue_max)
%         ctrs(j,1:2)=Blue_max(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
%     end
%     
%     for j=1:numel(Blue_max)
%         bbx(j,1:4)=Blue_max(j).BoundingBox;
%     end
%     for j=1:numel(image_data_r)
%         pxmean(j)=image_data_r(j).MeanIntensity; % loop goes thorugh each element and extracts centroid from structure B
%     end
%     for j=1:numel(Blue_max)
%         pxlist=Blue_max(j).PixelIdxList; % loop goes thorugh each element and extracts centroid from structure B
%         k=length(pxlist);
%         pxlist_z(1:k,j)=pxlist;
%     end
    %%
    %pxlist_z(pxlist_z==0)=nan();
%     Imax_z_filt=reshape(Imax_z,[],1);
%     Px_max_filt=reshape(Px_max,[],1);
%     size_max=size(pxlist_z);
%     for j=1:size_max(:,1)
%         for k=1:size_max(:,2)
%             pos=pxlist_z(j,k);
%             if pos~=0
%             Max_Ind(j,k)=Imax_z_filt(pos);
%             PX_max(j,k)=Px_max_filt(pos);
%             else Max_Ind(j,k)=nan();
%                 PX_max(j,k)=nan();
%             end
%         end
%     end
%     
%     Max_avg=mean(Max_Ind,1,'omitnan');
%     MI_avg=mean(PX_max,1,'omitnan');
%     cells_new=Max_avg*between_slice;
%     expresh_new=MI_avg;
%     cells_z_new=z_max-cells_new;
%     
%     cells_new_sh=reshape(cells_new,[],1);
%     expresh_new_sh=reshape(expresh_new,[],1);
%     s=numel(image_data);
%     
%     for k=1:s
%     cells_new_tot(k,position)=cells_new_sh(k,1);
%     expresh_new_tot(k,position)=expresh_new_sh(k,1);
%     end
%     cells_z_new_tot=z_max-cells_new_tot;
%     for j=1:numel(image_data_r)
%         Exp_r(j)=image_data_r(j).MeanIntensity; % identifies the particular pixel ID for each object. Identifies where the cells are
%     end
%     Red_new_sh=reshape(Exp_r,[],1);
%     
%     for k=1:numel(image_data_r)
%     Red_new_tot(k,position_r)=Red_new_sh(k,1);
%     end
%     %cells_z_new_tot=z_max-cells_new_tot;
%    
%    
% 
% else %cells_new_tot(1,position)=nan();
%      %expresh_new_tot(1,position)=nan();
%      %Green_new_tot(1,position_g)=nan();
%      Red_new_tot(1,position_r)=nan();
%  end
% 
%     Green_new_tot(Green_new_tot==0)=[];
%     Red_new_tot(Red_new_tot==0)=[];
%     cells_new_tot(cells_new_tot==0)=[];
%     Z_depth_data=reshape(cells_new_tot,[],1);
%    % expresh_new_tot(expresh_new_tot==0)=[];
%     Green_data=reshape(Green_new_tot,[],1);
%     Red_data=reshape(Red_new_tot,[],1);
% 
% %    end

%%
     
    %% gose throught each image in the z-stack
%     fprintf('Reading Individual Images from Image Stack %1.0f\n',position)
%     for i=1:num_z_stack(:,3)
%         file_s=strcat(dirname,image_list(position).name); % creates a variable with the full path to image # 1 (if there was a 2 inside the parenthesis, that would be image #2)
%         file=file_s{1};
%         info=imfinfo(file);
%         num_images = numel(info);
%         z_image=imread(file,i); % reads the image into variable im, witht he format it is saved in
%         
%         %% processing indavidual z-stack image
%         
%         z_rray_image=mat2gray(z_image,[0 2^16-1]); % normalizes the image from 0 to 1, the limits should be consistent with the bit depth of the image (2^12-1, in this case since 12 bit)
%         
%         z_filtered_image = bpass(z_gray_image,noise,object,thresh);
%         BW2z=imdilate(z_filtered_image,strel('sphere',2));
%         BW3z=bwareaopen(BW2z,100); %this is taking up a lot of time
%         BW4z=imfill(BW3z,'holes');
%         BW5z=imclearborder(BW4z);
%         
%         z_image_data_o=regionprops(BW5z,z_gray_image,'Centroid','BoundingBox','MeanIntensity'); %this is taking up a lot of time change
%        %%
%        unprocessed_imagez=z_image;
% Bz=imresize(unprocessed_imagez, 0.25);
% %imshow(B)
% Cz=rescale(Bz);%not BW image
% %figure(1),imshowpair(B,C,'montage')
% %im=C(:,:,:) %BW1 is black when used this
% %[L A Bl] = imsplit(rgb2lab(C))
% blueChannelz = Cz;%this is a BW picture
% %figure(1),imshow(blueChannel)
% %figure(2),histogram(blueChannel)
% blueChannel2z=blueChannelz>0.3; 
% %figure(3),imshow(blueChannel2) %why this is BW picture?
% BW1nz=mat2gray(blueChannel2z); %removed [0 2^8-1]) and it shows an area now %Q-can I show both L and b in Lab color space
% BW2nz=imdilate(BW1nz,strel('sphere',2)); % dilates objects
% BW3nz=imbinarize(BW2nz);
% %CC=bwconncomp(BW2);
% BW4nz = bwareafilt(BW3nz,[1 50]);
%     %figure(position);imshow(BW5,'InitialMagnification','fit')
%  %%   
%  z_image_data=regionprops(BW4nz,blueChannel2,'all');
%         
%         %% checking if cells are in this slice
%         if numel(z_image_data)~=0 && numel(image_data)~=0
%             
%             %%
%             for j=1:numel(z_image_data)
%                 ctrsz(j,1:2,i)= z_image_data(j).Centroid; % loop goes thorugh each element and extracts centroid from structure B
%             end
%             %% creats a 2D matrix cataloging the distance of centrods in this slice vs centrods in max image
%             clear distances
%             for j=1:numel(image_data)
%                 for k=1:numel(z_image_data)
%                     distances(j,k,1)=pdist([ctrs(j,:,1);ctrsz(k,:,i)]);
%                     
%                 end
%             end
%             
%             %% idenifies the min distance a centrod in the slice has to centrods in the max
%             distances_size=size(distances);
%             [min_distance(1:distances_size(:,2),i),I(1:distances_size(:,2),i)]=min(distances,[],1);
%             
%             
%             %% finds the area that the max and slice have in common
%             clear bbxz
%             clear area
%             for j=1:numel(z_image_data)
%                 bbxz(j,1:4)=z_image_data(j).BoundingBox;
%             end
%             area=rectint(bbx,bbxz);
%             
%             %% idenifies the mean intensity of each object in slice
%             clear pxmeanz
%             for j=1:numel(z_image_data)
%                 pxmeanz(j,i)=z_image_data(j).MeanIntensity; % loop goes thorugh each element and extracts centroid from structure B
%             end
%             %%  isoliites cells and indexs there position in matrix image
%             
%             %             clear Cells
%             %             area_size=size(area);
%             %             for j=1:area_size(:,1)
%             %                 for k=1:area_size(:,2)
%             %                     if area(j,k)~=0
%             %                          Cells(j,k)=distances(j,k);
%             %                          filtered_min_distance(j,i)=min_distance(j,i);
%             %                          pos(k,i)=k;
%             %                     end
%             %                 end
%             %             end
%             
%             clear Cells
%             image_size=size(area);
%             for j=1:image_size(:,1)
%                 for k=1:image_size(:,2)
%                     if area(j,k)~=0
%                         Cells(j,k)=min_distance(k,i);
%                     end
%                 end
%             end
%             %%
%             if exist('Cells','var')~=0
%                 %%
%                 image_sizez=size(Cells);
%                 for j=1:image_sizez(:,1)
%                     for k=1:image_sizez(:,2)
%                         if distances(j,k)==Cells(j,k)
%                             filtered_min_distance(j,i)=min_distance(k,i);
%                             filtered_index(j,i)=I(k,i);
%                             filtered_ctrs(j,i,1:2)=ctrsz(k,1:2,i);
%                             pos(k,i)=k;
%                         end
%                     end
%                 end
%                 %% conferms that cells were idendified in this slice
% %                 % maybe take out???
% %                                  if nnz(Cells)>length(pos)
% %                                      pos(k+1:nnz(Cells),i)=1;
% %                                  end
%                 %% isolates the mean pixle values for cells (odjects) in this slice
%                 filtered_min_size=size(filtered_min_distance);
%                 for j=1:filtered_min_size(:,1)
%                     if filtered_min_distance(j,i)~=0
%                         pxmean_filtered(j,i)=pxmean(j);
%                     end
%                 end
%                 
%                 %% puts the mean pixle values in a matrix for all slices
%                 counter_list(1:length(pxmean_filtered),1)=zeros();
%                 pxmean_filtered_size=size(pxmean_filtered);
%                 pos_size=size(pos);
%                 %%
%                 for j=1:pxmean_filtered_size(:,1)
%                     counter=1+nnz(counter_list);
%                     if counter<=numel(z_image_data) && counter<=pos_size(:,1)
%                         if pxmean_filtered(j,i)~=0 && pos(counter,i)~=0
%                             pxmean_filtered_z(j,i)=pxmeanz(counter,i);
%                             counter_list(j,1)=1;
%                         elseif pxmean_filtered(j,i)~=0
%                             counter_list(j,1)=1;
%                         end
%                         if numel(z_image_data)>counter
%                             counter=1+nnz(counter_list);
%                         end
%                     end
%                 end
%                 %% finds the Z depth where the a cell has its highest pixle value
%                 %and adds it to a final matrix with all cells z-depth
%                 if exist('pxmean_filtered_z')~=0
%                 [Max,Index]=max(pxmean_filtered_z,[],2);
%                 for k=1:length(Max)
%                     if Max(k,1)==0
%                         cells_z_depth(k,1)=nan();
%                     else cells_z_depth(k,1)=between_slice*Index(k,1);
%                     end
%                     cells_z_traveledn(k,1)=z_max-cells_z_depth(k,1);
%                     if cells_z_traveledn(k,1)>totaldepth
%                         cells_z_traveledn(k,1)=nan();
%                     end
%                     cells_z_traveled=rmmissing(cells_z_traveledn);
%                 end
%                 end
%             end
%         end
%     end
%         %%
%         if exist('cells_z_depth','var')~=0 && isempty(cells_z_traveled)~=1
%             clear ctrsa Ip Ip_1 Indexa
%             for j=1:numel(image_data)
%                 Indexa=Index(j,1);
%                 ctrsa(j,1:2,1)=filtered_ctrs(j,Indexa,1:2);
%                 ctrsa(j,3,1)=Indexa*between_slice;
%             end
%             
%             %%
%             disp(cells_z_traveled)
%             disp(size(cells_z_traveled))
%             %%
%             clear ctrsv dist
%             for j=1:numel(image_data)
%                 for k=1:numel(image_data)
%                 ctrsv=ctrsa(j,1:3);
%                 dist(j,k)=pdist([ctrsv(1,:);ctrsa(k,:)]);
%                 end
%             end
%             dist(dist==0)=NaN;
%             %%
%             avgdist(:,position)=mean(dist,2,'omitnan');
%             cells_z_traveledi(:,position)=cells_z_traveled(:,1);
%             %% writes the z-depth of each cell detected in to an excel file
%             cells_z_traveled1=reshape(cells_z_traveled,[],1);
%             avgdist1=reshape(avgdist,[],1);
%             row=length(cells_z_traveled);
%             str=['pos',int2str(position)];
%             heading={join(str)};
%             xlswrite(sheetname,heading(1),string(conditions(condition)),alpha(position)+"1:"+alpha(position)+"1")%writes the condition name
%             for r=1:row
%                 xlswrite(sheetname,cells_z_traveled(r),string(conditions(condition)),alpha(position)+""+(r+1)+":"+alpha(position)+""+(r+1))%writes the expression value
%             end
%             %adds data to final array that will contain all data for all positions
%             if position==1
%                 for i=1:length(cells_z_traveled)
%                     all_pos_data(i,1)=cells_z_traveled(i);
%                 end
%             length_all_pos(position,1)=length(all_pos_data);    
%             else
%                 for i=1:length(cells_z_traveled)
%                     all_pos_data(max(length_all_pos)+i,1)=cells_z_traveled(i);
%                 end
%             end
%             length_all_pos(position,1)=length(all_pos_data);
%             
%         end
%     end
%     deviation=std(all_pos_data);
%     mean_z=mean(all_pos_data);
%     cell_count=length(all_pos_data);
%     std_grror=deviation/sqrt(cell_count);
%     max_z=max(all_pos_data);
%     conditionname=string(conditions(condition));
%     xlswrite(sheetname,conditionname{1},'All Positions',alpha(condition)+"1:"+alpha(condition)+"1")%writes the condition name
%     %file=file_s{1};
%     xlswrite(sheetname,conditionname{1},'Summary',alpha(condition+1)+"1:"+alpha(condition+1)+"1")
%     disp('Writing all Positions')
%     for r=1:cell_count
%         xlswrite(sheetname,all_pos_data(r),'All Positions',alpha(condition)+""+(r+1)+":"+alpha(condition)+""+(r+1))%writes the expression value
%     end
%     xlswrite(sheetname,mean_z,'Summary',alpha(condition+1)+"2:"+alpha(condition+1)+"2")
%     xlswrite(sheetname,deviation,'Summary',alpha(condition+1)+"3:"+alpha(condition+1)+"3")
%     xlswrite(sheetname,cell_count,'Summary',alpha(condition+1)+"4:"+alpha(condition+1)+"4")
%     xlswrite(sheetname,std_grror,'Summary',alpha(condition+1)+"5:"+alpha(condition+1)+"5")
%     xlswrite(sheetname,max_z,'Summary',alpha(condition+1)+"6:"+alpha(condition+1)+"6")
%     for i=1:length(all_pos_data)
%     all_conditions(i,condition)=all_pos_data(i);
%     end
% end
% if(condition==2)
% [h,p] = ttest2(all_conditions(:,1),all_conditions(:,2));
% xlswrite(sheetname,stats(6),'Summary',"A7:A7")
% xlswrite(sheetname,p,'Summary',"B7:B7")
% end
% disp('done')
