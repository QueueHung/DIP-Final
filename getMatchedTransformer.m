function [L,a,b] = getMatchedTransformer( target_image , guided_image , target_weight , guided_weight , label , guided_ln_x_)
    % Some const
    SKY = 0;
    COLOR_L = 1;
    COLOR_a = 2;
    COLOR_b = 3;

    % L = sum( wn_x_ .* tn_x_ , 'all');
    target_ln_x_ = getMean(target_image,target_weight.predict_label,label,COLOR_L);
    % L =  target_image(:,:,1) + guided_ln_x_(label + 1) - target_ln_x_;
    L = guided_ln_x_(label + 1) - target_ln_x_;

    % Get chroma of target and guided 
    [ mean_target , cov_target ] = getChrominance(label,target_weight.predict_label,target_image);
    [ mean_guided , cov_guided ] = getChrominance(label,guided_weight.predict_label,guided_image);


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