% NOTE : if label is not in the image , weight must be ZERO
function rtv_weights = getNormalizeWeight(weights,threshold,isInTarget)% label : 因為天空的weight要為0

    [ lbl_deep , rows , cols] = size(weights);
    rtv_weights = weights;

    % disp(size(weights));

    for col = 1 : cols
        for row = 1 : rows 
            
            for i = 1 : lbl_deep
                if( rtv_weights(i,row,col) < threshold ) 
                    rtv_weights(i,row,col) = 0;
                end
                if( ~isInTarget(i) )
                    rtv_weights(i,row,col)= 0 ;
            end
            
            % rtv_weights(label+1,row,col) = 0;

            sm = sum(rtv_weights(:,row,col));
            if(sm == 0)
                sm = 1;
            end
            rtv_weights(:,row,col) = rtv_weights(:,row,col)/sm;

        end
    end


    % disp('Test');
    
end