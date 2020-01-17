function [matrixInv] = invreg(matrix, alpha)
    % Performs an inversion employing an SVD with Thikonov Regularization
    % on the singular values
    if nargin < 2 || isempty(alpha)
        alpha = 2;
    end
    
    [U,S,V] = svd(matrix); % U*S*V' = matrix
    S = diag(S);
    S = S.^2./(S.^2 + alpha)./S;
    S = diag(S);
    
    matrixInv = V*S*U';
end