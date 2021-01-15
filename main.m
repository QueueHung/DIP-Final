clc; clear all; close all;
% Some const
SKY = 0;
COLOR_L = 1;
COLOR_a = 2;
COLOR_b = 3;

target_name = '0137';
guided_name = '0005';

target_weight   = load(['MAToutput/',target_name,'.mat']);
guided_weight   = load(['MAToutput/',guided_name,'.mat']);

image_size = size(imread(['database/image/',target_name,'.png']));
target_image     = rgb2lab( imresize( sky_replace_v2(str2num(target_name), str2num(guided_name)) , [500 500] ) );
target_ori_image = rgb2lab( imresize( imread(['database/image/',target_name,'.png']) , [500 500] ) );
guided_image     = rgb2lab( imresize( imread(['database/image/',guided_name,'.png']) , [500 500] ) );

figure(str2double(target_name));
imshow(lab2rgb(target_ori_image));

figure(str2double(guided_name));
imshow(lab2rgb(guided_image));

threshold = 1/14;


% Apply Color Transfer function on Target Image

% transformFunc = @(x) getTransformedImage( x ,target_image,guided_image , target_weight,guided_weight , i,j , threshold);
% img_rlt = arrayfun( transformFunc , target_image);

% imshow(img_rlt);
isInTarget = false(1,14);
for i = 0:13
    if(  ismember(i, target_weight.predict_label ) )
        isInTarget(i+1) = true;
    end
end
isInGuided = false(1,14);
for i = 0:13
    if(  ismember(i, guided_weight.predict_label ) )
        isInGuided(i+1) = true;
    end
end


% Get mean for L of sky
mean_l_sky_target = getMean(target_image,target_weight.predict_label,SKY,COLOR_L);
mean_l_sky_guided = getMean(guided_image,guided_weight.predict_label,SKY,COLOR_L); 

% Get mean for a of sky
mean_a_sky_target = getMean(target_image,target_weight.predict_label,SKY,COLOR_a);
mean_a_sky_guided = getMean(guided_image,guided_weight.predict_label,SKY,COLOR_a); 

% Get mean for b of sky
mean_b_sky_target = getMean(target_image,target_weight.predict_label,SKY,COLOR_b);
mean_b_sky_guided = getMean(guided_image,guided_weight.predict_label,SKY,COLOR_b); 

% mean of color
color_diff = sqrt( ( mean_a_sky_target - mean_a_sky_guided )^2 + ( mean_b_sky_target - mean_b_sky_guided )^2 );

% Get beta
% TODO : 'Color' in paper is a and b
% TODO : Beta is mean of ( mean of a ) and ( mean of b )
cBeta = tanh(color_diff);

% Get Ln(x) : Mean of every transformed L
ln_x_ = zeros(1,14,'single'); 
for i = 0:13
    LIn = getMean(target_ori_image,target_weight.predict_label,i,COLOR_L);
    LRn = getMean(guided_image,guided_weight.predict_label,i,COLOR_L);
    ln_x_(i+1) = LIn + cBeta*( LRn - LIn );
end

% Get foreground's new L
target_img_l = target_ori_image(:,:,COLOR_L);
mean_without_sky_target_l = mean( target_img_l( target_weight.predict_label ~= SKY ) , 'all' );
target_img_a = target_ori_image(:,:,COLOR_a);
mean_without_sky_target_a = mean( target_img_a( target_weight.predict_label ~= SKY ) , 'all' );
target_img_b = target_ori_image(:,:,COLOR_b);
mean_without_sky_target_b = mean( target_img_b( target_weight.predict_label ~= SKY ) , 'all' );

guided_img_l = guided_image(:,:,COLOR_L);
mean_without_sky_guided_l = mean( guided_img_l( guided_weight.predict_label ~= SKY ) , 'all' );
guided_img_a = guided_image(:,:,COLOR_a);
mean_without_sky_guided_a = mean( guided_img_a( guided_weight.predict_label ~= SKY ) , 'all' );
guided_img_b = guided_image(:,:,COLOR_b);
mean_without_sky_guided_b = mean( guided_img_b( guided_weight.predict_label ~= SKY ) , 'all' );

color_diff_without_sky = sqrt( ( mean_without_sky_guided_a - mean_without_sky_target_a )^2 + ( mean_without_sky_guided_b - mean_without_sky_target_b )^2 );

nBeta = tanh(color_diff_without_sky);

L_non_matched = mean_without_sky_target_l + nBeta*( mean_without_sky_guided_l - mean_without_sky_target_l );

img_rlt = target_image;
%figure(100);
%imshow(lab2rgb(img_rlt));

wn_x_ = getNormalizeWeight(target_weight.predict_value,1/(sum(isInTarget) - 1), isInTarget);
% for lbl = 1:13
%     if(~ismember(lbl, target_weight.predict_label))
%         continue;
%     end
%     if(isInGuided(lbl + 1))
%         [ L , a , b ] = getMatchedTransformer(target_image , guided_image , target_weight , guided_weight , lbl , ln_x_);
%     else
%         [ L , a , b ] = getNOTMatchedTransformer(target_image , guided_image , target_weight , guided_weight , lbl , repmat(L_non_matched,[1 14]));
%     end

%     likelihood = wn_x_(lbl + 1,:,:);
%     likelihood = likelihood(1,:,:);
%     img_rlt(:,:,1) = img_rlt(:,:,1) + likelihood(1).*L;
%     img_rlt(:,:,2) = img_rlt(:,:,2) + likelihood(1).*a;
%     img_rlt(:,:,3) = img_rlt(:,:,3) + likelihood(1).*b;
%     % figure(lbl);
%     % imshow(lab2rgb(img_rlt));
% end
%target_img = rgb2lab(sky_replace_v2(str2num(target_name), str2num(guided_name)))
img_rlt = getBasicColorTransformer( target_image , guided_image , target_weight , guided_weight);

% for i = 1 : 500
%     for j = 1:500
%         % Get Wn_x
%         wn_x_ = getNormalizeWeight(target_weight.predict_value(:,i,j).',0.1,SKY);

%         color_target = target_image(i,j,:);
%         color_target = color_target(:).';
        
%         if(isInGuided( target_weight.predict_label(i,j) + 1 ))
%             % Use T(x)
%             img_rlt(i,j,:) = getMatchedTransformer( color_target , target_image , guided_image , target_weight , guided_weight , target_weight.predict_label(i,j) , tn_x_ , wn_x_ );
%         else
%             % Use Color Temperature transform
%             img_rlt(i,j,:) = getNOTMatchedTransformer( color_target , guided_image , guided_weight , L_non_matched ) ;
%         end 
%     end
% end

figure(1000);
img_rlt =  imresize(lab2rgb(img_rlt), image_size(1:2));
imshow( img_rlt );
imwrite( img_rlt, "ImageOutput/Output.jpg");

