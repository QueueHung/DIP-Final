function [clr_mean, cov_mat] = getChrominance(label, predict_label, img)
	%Returns chroma of an image segment as num_pixels x 2 matrix
	img_a = img(:,:,2);
	img_b = img(:,:,3);
	img_chroma_a = img_a((predict_label == label ));
	img_chroma_b = img_b((predict_label == label ));
    
    img_chroma = [img_chroma_a img_chroma_b];

	clr_mean    = mean(img_chroma);
	cov_mat     = cov(img_chroma);


end