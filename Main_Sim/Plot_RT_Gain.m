function fig = Plot_RT_Gain(gain, RTs)
    fig = figure;
    round_trips = linspace(1, RTs, RTs);
    plot(round_trips,gain,'b','LineWidth',1.5)
    title('Gain Over Propagation')
    xlabel('Gain medium pass'); ylabel('Gain (P_{after} / P_{before})');
end