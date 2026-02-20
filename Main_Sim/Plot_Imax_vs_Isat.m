function fig = Plot_Imax_vs_Isat(outputs, I_sat, Nhit, zs)
    fig = figure;
    s = [0, cumsum(abs(diff(zs)))];
    hit_idx = Nhit:Nhit:numel(zs);
    s_hits = s(hit_idx);
    plot(s_hits, outputs.Imax, 'DisplayName','Max intensity','LineWidth',1.5); hold on;
    yline(I_sat, '--', 'DisplayName','I_{sat}','LineWidth',1.5);
    xlabel('Total propagated distance [m]');
    ylabel('Intensity [W/m^2]');
    title('Max Intensity Over Propagation');
    legend('Location','best');
    grid on;
end