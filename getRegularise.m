function rtv_mat = getRegularise(mat,reg)
    rtv_mat = mat;
    if rtv_mat(1,1)<reg
        rtv_mat(1,1)=reg;
    end
    
    if rtv_mat(2,2)<reg
        rtv_mat(2,2)=reg;
    end

end