function [loss1, loss2] = Analyze_Loss(laser, mirror)
    loss1 = exp(-2*(mirror(1).D/2)^2 / laser.w0); % loss from finite aperture
    loss2 = exp(-2*(mirror(2).D/2)^2 / laser.w0);
end