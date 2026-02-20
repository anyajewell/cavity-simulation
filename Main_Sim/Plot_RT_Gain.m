function fig = Plot_RT_Gain(gain)

    fig = figure;

    n = numel(gain);
    passes = 1:n;

    plot(passes, gain, 'b', 'LineWidth',1.5)

    title('Gain Over Propagation')
    xlabel('Gain medium pass')
    ylabel('Gain (P_{after} / P_{before})')
    grid on;

end