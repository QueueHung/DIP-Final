function [L,a,b] = getNOTMatchedTransformer( target_image , guided_image , target_weight , guided_weight , label , guided_ln_x_ )
    % Some const
    SKY = 0;
    COLOR_L = 1;
    COLOR_a = 2;
    COLOR_b = 3;

    % L = sum( wn_x_ .* tn_x_ , 'all');
    img_l = target_image(:,:,1);
    target_ln_x_ = mean( img_l( target_weight.predict_label ~= SKY ) , 'all' );
    % L =  target_image(:,:,1) + guided_ln_x_(label + 1) - target_ln_x_;
    L = guided_ln_x_(label + 1) - target_ln_x_;

    % Get chroma of target and guided 
    [ mean_target , cov_target ] = getChrominance(label,target_weight.predict_label,target_image);

    % Guided chroma
    % Returns chroma of an image segment as num_pixels x 2 matrix
	img_a = guided_image(:,:,2);
	img_b = guided_image(:,:,3);
	img_chroma_a = img_a((guided_weight.predict_label ~= SKY ));
	img_chroma_b = img_b((guided_weight.predict_label ~= SKY ));
    
    img_chroma = [img_chroma_a img_chroma_b];

	mean_guided    = mean(img_chroma);
    cov_guided     = cov(img_chroma);
    
    % [ mean_guided , cov_guided ] = getChrominance(label,guidde_weight.predict_label,guided_image);


    % Regularise CovMat 
    cov_target = getRegularise(cov_target,7.5);
    cov_guided = getRegularise(cov_guided,7.5);

    % Calculate T
    T = cov_target^-0.5*( cov_target^0.5 * cov_guided * cov_target^0.5 )^0.5;

    %
    x1 = T(1,1);
    x2 = T(1,2);
    x3 = T(2,1);
    x4 = T(2,2);
    a = target_image(:,:,2) - mean_target(1);
    b = target_image(:,:,3) - mean_target(2);
    a = (a*x1) + (b*x2) + mean_guided(1);
    b = (a*x3) + (b*x4) + mean_guided(2);

    a = a - target_image(:,:,COLOR_a);
    b = b - target_image(:,:,COLOR_b);
    % [a b] = T*( color_target(2:3) - [ mean_target_a mean_target_b ] ) + [mean_guided_a mean_guided_b ];


end