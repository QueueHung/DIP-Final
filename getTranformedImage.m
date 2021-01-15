function rtv_img = getTransformedImage( x , target_image,guided_image , target_weight,guided_weight , i,j , threshold  )
    
    
    if(ismember( x , guided_weight.predict_label))
        % Use T(x)
        rtv_img = getMatchedTransdformer(target_image,guided_image , target_weight,guided_weight , i,j , threshold );
    else
        % Use Color Temperature transform
        rtv_img = getNOTMatchedTransdformer(target_image,guided_image , i,j );
    end 

end