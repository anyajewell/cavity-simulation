function [p_out,q_out,d_out,d_us_out] = output(x_in,s_in,loslength)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% x_in = states [theta_s; omega_s; theta_g_rel; omega_g_rel] 
% s_in = LOS bearing angle trajectory of the target relative to the pursuer (rad) 
% los_length = distance between target/pursuer satellites (m)
%
% OUTPUTS:
% p_out = quantized target satellite angle relative to pursuer (quantized s_in)
% q_out = quanitzed pursuer satellite angle (quantized x_in)
% d_out = fine pointing error
% d_us_out = fine pointing error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Quantization Parameters
% Set quantization for coarse sensing of q
% set denom of quantizionq to be 1.
quantizationq = deg2rad(0.01)/1024;  % Set quantization = input range/desired quantization
% Set quantization for coarse sensing of p
quantizationp = deg2rad(0.5)/1024;  % Set quantization = input range/desired quantization

%% P_out
% p_out     = quantizationp .* round( (s_in(1,:)./loslength )./quantizationp);  % This accomplishes the quantization  (if skc is in m)
p_out     = quantizationp .* round( (s_in(1,:) )./quantizationp);  % This accomplishes the quantization (if skc is in rad) (if skc is in rad I need to remove loslength)

%% q_out
q_out     = quantizationq .* round( (x_in(1,:) )./quantizationq); % This accomplishes the quantization 

%% d_us_out
% Adding together x_in(1,:) + x_in(3,:) 
d_us_out  = s_in(1,:) - ( x_in(1,:) + x_in(3,:) ); % p - q but clean and including gimbal angle (if skc is in rad)
% d_us_out  = s_in(1,:)./loslength - ( x_in(1,:) + x_in(3,:) ); % p - q but clean and including gimbal angle (if skc is in m)

%unsaturated version
d_out = d_us_out;

% % Saturate d measurements to model sensor limitations
% if abs(d_us_out) > deg2rad(0.5)
%     d_out = sign(d_us_out)*deg2rad(0.5);
% else
%     d_out = d_us_out;
% end

end