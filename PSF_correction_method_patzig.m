% PSF correction method
function [reconstructed_object] = PSF_correction_method_patzig(image, fieldMap, t2starMap, params, scan)

reconstructed_object = zeros(size(image));

parfor ic1=1:size(image,2)
    
    PSFmatrix_parfor    = zeros(size(image));
    fieldMapCol_parfor  = fieldMap(:,ic1);
    t2starMapCol_parfor = t2starMap(:,ic1);
    for ic2=1:length(fieldMapCol_parfor)
        PSF = getPSF(params, fieldMapCol_parfor(ic2), t2starMapCol_parfor(ic2), scan);
        PSFmatrix_parfor(:,ic2) = circshift(PSF,-ceil(size(image,2)/2)+ic2);
    end
    if isreal(image)
        PSFmatrix_parfor = abs(PSFmatrix_parfor);
        PSFmatrix_parfor = PSFmatrix_parfor./repmat(abs(sum(PSFmatrix_parfor)), size(image,1), 1);
    end
    
    inv_matrix = invreg(PSFmatrix_parfor, params.alpha);
    reconstructed_object(:,ic1) = inv_matrix*image(:,ic1);
end

end
