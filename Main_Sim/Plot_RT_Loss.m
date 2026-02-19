function fig = Plot_RT_Loss(loss_frac, RTs)
    fig = figure;
    round_trips = linspace(1, RTs, RTs);
    plot(round_trips,loss_frac,'b','LineWidth',1.5)
    title('Loss Fraction Over Propagation')
    xlabel('Round trip number'); ylabel('Loss fraction');
end