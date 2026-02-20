% function fig = Plot_Mirror_Loss(loss1, loss2, Nhit, zs)
%     fig = figure;
%     s = [0, cumsum(abs(diff(zs)))]; % total distance traveled at each step (same length as zs)
%     hit_idx = Nhit:Nhit:numel(zs); % steps where mirror hits occur
%     m1_idx = hit_idx(1:2:end);
%     m2_idx = hit_idx(2:2:end);
%     s_m1 = s(m1_idx);
%     s_m2 = s(m2_idx);
%     plot(s_m1, loss1, 'b.','MarkerSize', 15,'DisplayName','Mirror 1 (RHS)'); hold on;
%     plot(s_m2, loss2, 'r.','MarkerSize', 15,'DisplayName','Mirror 2 (LHS)');
%     title('Loss Fraction Over Total Propagated Distance');
%     xlabel('Total propagated distance [m]');
%     ylabel('Loss fraction');
%     legend('Location','best');
% end

function fig = Plot_Mirror_Loss(loss1, loss2, zs, L)

    fig = figure;

    % Total propagated distance
    s = [0, cumsum(abs(diff(zs)))];

    % Tolerance for floating point comparison
    tol = 1e-9 * max(1, abs(L));

    % Find mirror hit indices from zs
    m1_idx = find(abs(zs - L/2) < tol);   % RHS mirror
    m2_idx = find(abs(zs + L/2) < tol);   % LHS mirror

    % Convert to propagated distance
    s_m1 = s(m1_idx);
    s_m2 = s(m2_idx);

    % Match lengths safely
    n1 = min(length(s_m1), length(loss1));
    n2 = min(length(s_m2), length(loss2));

    plot(s_m1(1:n1), loss1(1:n1), 'b.', ...
         'MarkerSize',15,'DisplayName','Mirror 1 (RHS)');
    hold on;

    plot(s_m2(1:n2), loss2(1:n2), 'r.', ...
         'MarkerSize',15,'DisplayName','Mirror 2 (LHS)');

    title('Loss Fraction Over Total Propagated Distance');
    xlabel('Total propagated distance [m]');
    ylabel('Loss fraction');
    legend('Location','best');
    grid on;

end