function fig = Plot_Loss_Over_Prop(loss_frac)
    fig = figure;
    n_mirror_hits = linspace(1, size(loss_frac));
    title({'Loss Fraction in x Over Propagation'})
    plot(n_mirror_hits,loss_fraction,'b','LineWidth',1.5)
    xlabel('Mirror hit number [m]'); ylabel('Loss fraction [m]');
    legend();
end