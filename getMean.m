function rtv_means = getMean(img,predict_label,label,color_space)
    img_l = img(:,:,color_space);
    rtv_means = mean( img_l( predict_label == label ) , 'all' );
    if(isnan(rtv_means))
        rtv_means = 0;
    end
end