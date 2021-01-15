
%addpath(genpath('custom_toolboxes'))
%atom_path = '/Neutron9/anjali.shenoy/dip_project/';
I = dir('dataset/image/*.png');
M = dir('dataset/mask/*.png');

no_of_images = size(I,1);
%w = waitbar(0,'Computing max sky areas')
for k = 1:no_of_images
    if k==144
        continue;
    end
    source = imread(['dataset/image/',M(k).name]);
    source_sky_mask = im2bw(imread(['dataset/mask/',M(k).name]));

    if isempty(find(source_sky_mask==0))
        max_rects(k).max_source_region = 0;
       continue;
    end
    [S,~,~,max_source_region,sr1,sr2,sc1,sc2] = FindLargestRectangles((source_sky_mask),source);
    max_rects(k).max_source_region = max_source_region;
    max_rects(k).index = [sr1,sr2,sc1,sc2];
    max_rects(k).region = S;
    k
    %waitbar(k/no_of_images,sprintf('percentage = %2.2f',(k*100)/no_of_images))
end
save('mat_files/max_rects.mat','max_rects')

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