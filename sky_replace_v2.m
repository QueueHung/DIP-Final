function I_replaced = sky_replace_v2(file1, file2)
    [input_region,I_replaced,reference_Im,reference_region,norm_F,ref_im_name] = replace_sky(file1, file2);


    %show
    I = dir('database/image/*.png');
    input = imread(['database/image/',I(file1).name]);
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
end
function[rough_mask,I_rep,I_ref,ref_region,norm_F,ref_im_name] = replace_sky(file1, file2)
%addpath(genpath('custom_toolboxes'))
    I = dir('database/image/*.png');
    P = dir('database/mask/*.png');  
    L = dir('MAToutput/*.mat');

    A = imread(['database/image/',I(file1).name]);
    target_sky_mask = imread(['database/mask/',P(file1).name]);
    target_sky_mask = im2bw(target_sky_mask);
    target = imresize(A, [500 500]);
    load (['MAToutput/', L(file1).name]);
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
    %%%%%%%%%% choose up to 5 skies to display %%%%%%%%%%%
    %w = waitbar(0,'Finding candidate skies')
    [~,~,~,max_target_region,tr1,tr2,tc1,tc2] = FindLargestRectangles((target_sky_mask),A);
    source = imread(['database/image/',I(file2).name]);
    source_sky_mask = im2bw(imread(['database/mask/',P(file2).name]));
    ref_im_name = I(file2).name;
    %%%%%%%%% replace sky using the chosen image %%%%%%%%%%%
    [S,~,~,max_source_region,sr1,sr2,sc1,sc2] = FindLargestRectangles((source_sky_mask),source);
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
    load (['MAToutput/', L(file2).name]);
    ref_region = predict_label;
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
