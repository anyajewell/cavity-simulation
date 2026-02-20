function fig = Plot_Imax_vs_Isat(outputs, I_sat, sim)
    fig = figure;
    s = [0, cumsum(abs(diff(outputs.zs)))];
    [hit_idx, ~, ~] = Get_Mirror_Hit_Idx(outputs.zs, sim.L);  % indices where either mirror hit occurred
    s_hits = s(hit_idx);
    plot(s_hits, outputs.Imax, 'DisplayName','Max intensity','LineWidth',1.5); hold on;
    yline(I_sat, '--', 'DisplayName','I_{sat}','LineWidth',1.5);
    xlabel('Total propagated distance [m]');
    ylabel('Intensity [W/m^2]');
    title('Max Intensity Over Propagation');
    legend('Location','best');
    grid on;
end