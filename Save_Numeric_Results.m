function [Es, wzs] = Save_Numeric_Results(E, Es, wzs, numSteps, H, X, Y, dx, dy)    

    for n = 1:numSteps
        
        if n == 1 % save values at the mirror face before propagating
            % no propagation yet
        else
            % propagation    
            FE = fft2(E); % transform beam to frequency domain
            FE = FE.*fftshift(H); % propagate beam in frequency domain 
            E = ifft2(FE); % transform back to space domain 
            
        end

        % Save variables of interest
        Es{n} = E; % save intermediate field
        I_num = abs(E).^2; % numeric intensity
        Pnum = sum(I_num(:)) * dx * dy; % total power
        xbar = sum(sum(I_num .* X)) * dx * dy / Pnum; % first moment x (x centroid)
        ybar = sum(sum(I_num .* Y)) * dx * dy / Pnum; % second moment y (y centroid)
        sigma_r_num = sqrt(sum(sum(I_num .* ((X-xbar).^2 + (Y-ybar).^2))) * dx * dy / Pnum); % second moment (variance)
        wzs(n) = 2 * sigma_r_num; % 1/e^2 waist at this point

        % Peak intensity
        %I0s(n) = max(I_num(:)); 

    end

end
