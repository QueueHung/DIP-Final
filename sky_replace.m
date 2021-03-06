

file1 = 4;
[input_region,I_replaced,reference_Im,reference_region,norm_F,ref_im_name] = replace_sky(file1);


%show
I = dir('dataset/image/*.png');
input = imread(['dataset/image/',I(file1).name]);
figure;
subplot(1,3,1)
imshow(input)
title('input')
subplot(1,3,2)
imshow(reference_Im)
title('reference')
subplot(1,3,3)
imshow(I_replaced)
title('output')
function[rough_mask,I_rep,I_ref,ref_region,norm_F,ref_im_name] = replace_sky(file1)
     %%%% MAT FILES %%%%
    load(['mat_files/','descriptor_212']); 
    load(['mat_files/','max_rects']); 
%addpath(genpath('custom_toolboxes'))
    I = dir('dataset/image/*.png');
    P = dir('dataset/mask/*.png');  
    L = dir('fcn_dat/*.mat');

    A = imread(['dataset/image/',I(file1).name]);
    target_sky_mask = imread(['dataset/mask/',P(file1).name]);
    target_sky_mask = im2bw(target_sky_mask);
    im_index = 1; %check filename
    threshold = 1.9;
    no_of_images = size(descriptor,2);
    target = imresize(A, [500 500]);
    %[~,rough_mask,res] = scene_parse(target);
    load (['fcn_dat/', L(file1).name]);
    rough_mask = predict_label;
    F = predict_value;
    for i=1:size(target,1)
        for j=1:size(target,2)
            m1 = min(F(:,i,j));
            m2 = max(F(:,i,j));
            norm_F(:,i,j) = (F(:,i,j) - m1) / (m2 - m1) ;
            sum1 = sum(norm_F(:,i,j));
            norm_F(:,i,j) = norm_F(:,i,j)/sum1;
        end
    end

    %%%%%%%%%% Compute target's descriptor '%%%%%%%%%%%
    H = [];
    sz = size(target);
    m = prod(sz(1:2));
    for gx = 1:3
        for gy = 1:3
            x = (floor(sz(1)/3))*gx;
            start_x = ((floor(sz(1)/3))*(gx-1)) + 1;
            y = (floor(sz(2)/3))*gy;
            start_y = ((floor(sz(2)/3))*(gy-1)) + 1;
            grid = norm_F(:,start_x:x,start_y:y); %check filename
            h = zeros(1,14);
            for label = 1:14
                h(label) = (1/m)*sum(sum(grid(label,:,:)));
            end
            H = [H; h];
        end
    end
    %global_histFig = figure('Name','Histogram','NumberTitle','off','Visible','off');
    %global_hist = histogram(target(:,:,:),14,'BinLimits',[0 1],'Normalization','probability','Visible','off');
    global_hist = histcounts(target(:,:,:),14);
    global_hist = global_hist/max(max(global_hist));
    H = [H; global_hist];
    target_descriptor = H;

    %%%%%%%%%% choose up to 5 skies to display %%%%%%%%%%%
    %w = waitbar(0,'Finding candidate skies')
    chosen_sky = 1;
    no_of_chosen_skies = 5;
    visited = zeros([1,no_of_images]);
    jk = 1; end_flag = 0;
    options = [];
    [~,~,~,max_target_region,tr1,tr2,tc1,tc2] = FindLargestRectangles((target_sky_mask),A);
    t_width = tc2-tc1;
    t_height = tr2-tr1;
    P_ta = t_width/t_height;
    P_ts = t_width*t_height;
    end_flag = 0; chosen_sky = 1;
    visited = zeros([1,no_of_images]);
    semantic_similarity = [];
    found = 0;
    for k = 1:no_of_images
        if k==144
            continue;
        end
        val = norm(descriptor(k).desc - target_descriptor);
        semantic_similarity = [semantic_similarity ; val];
    end

    [ASorted AIdx] = sort(semantic_similarity);

    for jk = 2:no_of_images
        k = AIdx(jk);
        source = imread(['dataset/image/',P(k).name]);
        source_sky_mask = im2bw(imread(['dataset/mask/',P(k).name]));

        if (max_rects(k).max_source_region == 0)
           continue;
        end
        S = max_rects(k).region; 
        max_source_region = max_rects(k).max_source_region;
        sr1 = max_rects(k).index(1);
        sr2 = max_rects(k).index(2);
        sc1 = max_rects(k).index(3); 
        sc2 =  max_rects(k).index(4);
        %% Compute Aspect Ratio and Resolution
        s_width = sc2-sc1;
        s_height = sr2-sr1;
        P_sa = s_width/s_height;
        P_ss = s_width*s_height;
        Q_s = min(P_ts,P_ss)/max(P_ts,P_ss);
        Q_a = min(P_ta,P_sa)/max(P_ta,P_sa);
        if Q_s>0.5 && Q_a>0.5
            options(chosen_sky)=k;
            chosen_sky = chosen_sky + 1;
        end
        if chosen_sky >= 5
            break
        end
    end
    if chosen_sky~=1
        [index]=select_image(4,options,I);
    end
    source = imread(['dataset/image/',I(options(index)).name]);
    ref_im_name = I(options(index)).name;
    %%%%%%%%% replace sky using the chosen image %%%%%%%%%%%
    S = max_rects(options(index)).region;
    [T,tmax_area,tr1,tr2,tc1,tc2]  = find_largest_convex_hull(target_sky_mask,A);
    source_sky = imresize(S,[size(T,1) size(T,2)]);
    I_rep = A;
    for i = 1:size(I_rep,1)
        for j=1:size(I_rep,2)
            if(target_sky_mask(i,j)==1)
                I_rep(i,j,:) = source_sky(i,j,:);
            end
        end
    end
    I_ref = source;
    load (['fcn_dat/', L(options(index)).name]);
    ref_region = predict_label;
    %[~,ref_region,~] = scene_parse(I_ref);
    %figure;imshowpair(target,I_rep,'montage');str = sprintf('Original & Final Image');title(str);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
function[S,max_area,r1,r2,c1,c2]  = find_largest_convex_hull(mask,image)
	[rows,cols,~] = size(mask);
	
	end_flag = 0; col = 1; start = 0; last = 0;
	while end_flag~=1 && col<=cols
		column = mask(:,col);
		if start==0&&numel(find(column==1))~=0
			start = col;
		elseif start~=0&&numel(find(column==1))==0
			last = col-1;
			end_flag = 1;
		end
		col=col+1;
	end
	c2 = start;
	if last == 0
		last = cols;
	end
	c1 = last;
	
	end_flag = 0; start = 0; last = 0; row = 1;
	while row<=rows&&end_flag~=1
		row_vec = mask(row,:);
		if start==0&&numel(find(row_vec==1))~=0
			start = row;
		elseif start~=0&&numel(find(row_vec==1))==0
			last = row-1;
			end_flag = 1;
		end
		row=row+1;
	end
	r2 = start;
	r1 = last;

	S = image(r2:r1,c2:c1,:);
	max_area = prod(size(S,1),size(S,2));
end
%%
function[index]=select_image(k,options,I)
  hfig=figure;
  for i=1:k
    subplot(2,2,i);
    im = imread(['dataset/image/',I(options(i)).name]);
    h{i}.h = imshow(im);
    set(h{i}.h, 'buttondownfcn', {@loads_of_stuff,i});
  end

  function loads_of_stuff(src,eventdata,x)
    if get(src,'UserData')
        set(src,'UserData',0)
        title('');
    else
        set(src,'UserData',1)
        title('Selected');
        %[filename,user_canc]=imsave(src);
        %image=imread(filename);
    end
    fprintf('%s:\n',num2str(x));
    index = x;
  end
uiwait(hfig);
end
%%
function [C, H, W, M,r1,r2,c1,c2] = FindLargestRectangles(I, image, crit, minSize)
    if (nargin<3)
      crit = [1 1 0];
    end
    if (nargin<4)
      minSize = [1 1];
    end
    p = crit;
    [nR nC] = size(I);
    if (minSize(1)<1), minSize(1)= floor(minSize(1)*nR); end
    if (minSize(2)<1), minSize(2)= floor(minSize(2)*nC); end
    if (max(I(:)) - min(I(:))==1),
      S = FindLargestSquares(I);
    else
      S = I;
    end
    n = max(S(:));
    W = S; % make a carbon copy of the matrix data
    H = S;
    C = ((p(1)+p(2)) + p(3)*S) .* S; % p(1)*width + p(2)*height + p(3)*height*width for height=width=S;
    d = round((3*n)/4);
    minH = max(min(minSize(1), d),1);
    minW = max(min(minSize(2), d),1);
    hight2width = zeros(n+1,1);  % Store array with largest widths aviable for a given height
    for r = 1 : nR               % each row is processed independently
      hight2width(:) = 0;        % reset the List
      for c = nC: -1 : 1         % go through all pixels in a row right to left
        s = S(r,c);              % s is a size of a square with its corner at (r,c)
        if (s>0)                 % if pixel I(r,c) is true
          MaxCrit = C(r,c);      % initialize the Max Criteria using square
          for hight = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
            width = hight2width(hight); % look up width for a given hight
            width = max(width+1,s);
            hight2width(hight) = width;
            Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
            if (Crit>MaxCrit),   % check if it produces larger Criteria
              MaxCrit = Crit;    % if it does than save the results
              W(r,c)  = width;
              H(r,c)  = hight;
            end % if Crit
          end % for hight
          C(r,c)  = MaxCrit;
        end % if s
        hight2width((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
      end % for c
    end
    clear hight2width
    width2hight = zeros(n+1,1);  % Store array with largest widths aviable for a given height
    for c = 1 : nC               % each column is processed independently
      width2hight(:) = 0;        % reset the List
      for r = nR: -1 : 1         % go through all pixels in a column bottom to top
        s = S(r,c);              % s is a size of a square with its corner at (r,c)
        if (s>0)                 % if pixel I(r,c) is true
          MaxCrit = C(r,c);      % initialize the Max Criteria using square
          for width = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
            hight = width2hight(width); % look up hight for a given width
            hight = max(hight+1,s);
            width2hight(width) = hight;
            Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
            if (Crit>MaxCrit),   % check if it produces larger Criteria
              MaxCrit = Crit;    % if it does than save the results
              W(r,c)  = width;
              H(r,c)  = hight;
            end % if Crit
          end % for width
          C(r,c)  = MaxCrit;
        end % if s
        width2hight((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
      end % for r
    end
    CC = C;
    CC( H<minH | W<minW ) = 0; % first try to obey size restrictions
    [~, pos] = max(CC(:));
    if (isempty(pos)), [~, pos] = max(C(:)); end % but when it fails than drop them
    [r c] = ind2sub(size(C), pos);
    M = false(size(C));
    r1 = r;
    r2 = (r+H(r,c)-1);
    c1 = c;
    c2 = (c+W(r,c)-1);

    M( r1:r2, c1:c2 ) = 1;
    C = image( r1:r2, c1:c2, :);
end
%%
function S = FindLargestSquares(I)
    [nr nc] = size(I);
    S = double(I>0);
    for r=(nr-1):-1:1
      for c=(nc-1):-1:1
        if (S(r,c))
          a = S(r  ,c+1);
          b = S(r+1,c  );
          d = S(r+1,c+1);
          S(r,c) = min([a b d]) + 1;
        end
      end  
    end  
end
