function [F, G] = makeFG(theta_, nhat_, mugap_, m_, p_)
%makeFG creates the F and G matrices for an ARX model from theta vector
%coefficients. m and p are the number of inputs and outputs, and nhat is 
%the model order. F is a p x p x nhat matrix and G is a p x m x nhat+1
%matrix. The 3rd dimension represents the index starting from nhat, decreasing to 0
%(G is offset by 1 since its index ends at 0).
%Example
%Equation: 
% y_k = -F_{nhat-mugap}*y_{k-nhat+mugap} ... -F_{nhat}*y_{k-nhat} + G_{0}*u_{k} + ... + G_{nhat}*u_{k-nhat}
%Code:
% y_k = -F(:,:,nhat-mugap)*y_{k-nhat+mugap} ... -F(:,:,nhat)*y_{k-nhat} + G(:,:,1)*u_{k} + ... + G_(:,:,nhat+1)*u_{k-nhat}

%Split theta into F and G parts
thFlength = ( (nhat_-mugap_)*p_ ); 
thF_temp  = theta_(:,1:thFlength);
% Pad thF_temp with zeros where mu-gap is 
thF       = [ zeros(p_,p_*mugap_) thF_temp ];
thG       = theta_(:,thFlength+1:end);

%Reshape ID theta into F and G
G = reshape(thG, p_, m_, nhat_+1);
F = reshape(thF, p_, p_, nhat_); 

end


