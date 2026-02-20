% function fig = Plot_R(R1, R2, Nhit, zs)
%     fig = figure;
%     s = [0, cumsum(abs(diff(zs)))]; % total distance traveled at each step (same length as zs)
%     hit_idx = Nhit:Nhit:numel(zs); % steps where mirror hits occur
%     m1_idx  = hit_idx(1:2:end);
%     m2_idx  = hit_idx(2:2:end);
%     s_m1 = s(m1_idx);
%     s_m2 = s(m2_idx);
%     plot(s_m1, R1, 'b.','MarkerSize', 15,'DisplayName','Mirror 1 (RHS)'); hold on;
%     plot(s_m2, R2, 'r.','MarkerSize', 15,'DisplayName','Mirror 2 (LHS)');
%     title('Reflection Coefficient Over Total Propagated Distance');
%     xlabel('Total propagated distance [m]');
%     ylabel('Reflection coefficient');
%     legend('Location','best');
% end

function fig = Plot_R(R1, R2, zs, L)

    fig = figure;

    % Total propagated distance at each logged sample
    s = [0, cumsum(abs(diff(zs)))];

    % Tolerance for floating point comparisons
    tol = 1e-9 * max(1, abs(L));

    % Indices where the beam is at each mirror location
    m1_idx = find(abs(zs - L/2) < tol);   % mirror 1 (RHS)
    m2_idx = find(abs(zs + L/2) < tol);   % mirror 2 (LHS)

    % Map hits to propagated distance
    s_m1 = s(m1_idx);
    s_m2 = s(m2_idx);

    % Match array lengths safely
    n1 = min(numel(s_m1), numel(R1));
    n2 = min(numel(s_m2), numel(R2));

    plot(s_m1(1:n1), R1(1:n1), 'b.', 'MarkerSize', 15, ...
        'DisplayName','Mirror 1 (RHS)'); 
    hold on;

    plot(s_m2(1:n2), R2(1:n2), 'r.', 'MarkerSize', 15, ...
        'DisplayName','Mirror 2 (LHS)');

    title('Reflection Coefficient Over Total Propagated Distance');
    xlabel('Total propagated distance [m]');
    ylabel('Reflection coefficient');
    legend('Location','best');
    grid on;

end