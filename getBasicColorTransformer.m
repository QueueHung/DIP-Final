function rtv_img = getBasicColorTransformer(target_img , guided_img , target_weight , guided_weight)
    target_img_l = target_img(:,:,1);
    target_img_a = target_img(:,:,2);
    target_img_b = target_img(:,:,3);
    
    target_std_sky_l = std(target_img_l(target_weight.predict_label == 0));
    target_std_sky_a = std(target_img_a(target_weight.predict_label == 0));
    target_std_sky_b = std(target_img_b(target_weight.predict_label == 0));

    target_std_sky = [target_std_sky_l target_std_sky_a target_std_sky_b];

    guided_img_l = guided_img(:,:,1);
    guided_img_a = guided_img(:,:,2);
    guided_img_b = guided_img(:,:,3);

    guided_std_sky_l = std(guided_img_l(guided_weight.predict_label == 0));
    guided_std_sky_a = std(guided_img_a(guided_weight.predict_label == 0));
    guided_std_sky_b = std(guided_img_b(guided_weight.predict_label == 0));

    guided_std_sky = [guided_std_sky_l guided_std_sky_a guided_std_sky_b ];

    disp(target_std_sky);
    % disp(guided_std_sky);

    target_mean_sky_l = getMean(target_img,target_weight.predict_label,0,1);
    target_mean_sky_a = getMean(target_img,target_weight.predict_label,0,2);
    target_mean_sky_b = getMean(target_img,target_weight.predict_label,0,3);

    guided_mean_sky_l = getMean(guided_img,guided_weight.predict_label,0,1);
    guided_mean_sky_a = getMean(guided_img,guided_weight.predict_label,0,2);
    guided_mean_sky_b = getMean(guided_img,guided_weight.predict_label,0,3);

    delta_mean_l = target_img(:,:,1) - target_mean_sky_l;
    delta_mean_a = target_img(:,:,2) - target_mean_sky_a;
    delta_mean_b = target_img(:,:,3) - target_mean_sky_b;

    rtv_img(:,:,1) = ( guided_std_sky(1) / target_std_sky(1) )*delta_mean_l + target_mean_sky_l;
    rtv_img(:,:,2) = ( guided_std_sky(2) / target_std_sky(2) )*delta_mean_a + target_mean_sky_a;
    rtv_img(:,:,3) = ( guided_std_sky(3) / target_std_sky(3) )*delta_mean_b + target_mean_sky_b;


end