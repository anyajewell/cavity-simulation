function fig = Plot_Mirror_Loss(loss1, loss2, Nhit, zs)
    fig = figure;
    s = [0, cumsum(abs(diff(zs)))]; % total distance traveled at each step (same length as zs)
    hit_idx = Nhit:Nhit:numel(zs); % steps where mirror hits occur
    m1_idx  = hit_idx(1:2:end);
    m2_idx  = hit_idx(2:2:end);
    s_m1 = s(m1_idx);
    s_m2 = s(m2_idx);
    plot(s_m1, loss1, 'b.','MarkerSize', 15,'DisplayName','Mirror 1 (RHS)'); hold on;
    plot(s_m2, loss2, 'r.','MarkerSize', 15,'DisplayName','Mirror 2 (LHS)');
    title('Loss Fraction Over Total Propagated Distance');
    xlabel('Total propagated distance [m]');
    ylabel('Loss fraction');
    legend('Location','best');
end