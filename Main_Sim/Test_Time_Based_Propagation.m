function tests = Test_Time_Based_Propagation
    tests = functiontests(localfunctions);
end

function testElapsedTimeAdvancesAcrossBothDirections(testCase)
    consts.c = 3e8;
    consts.eps0 = 8.8541878128e-12;

    sim.N = 16;
    sim.Lx = 1;
    sim.dx = sim.Lx/sim.N;
    sim.dt = 0.4e-9;
    sim.direction = 1;
    [sim.x, sim.y, sim.X, sim.Y] = Create_Grid(sim.N, sim.Lx, sim.dx);

    laser.Ld = 1e-6;
    laser.k0 = 2*pi/laser.Ld;
    laser.pos = 0;
    laser.time = 0;
    laser.Gau = exp(-(sim.X.^2 + sim.Y.^2)/0.1);
    laser.P_ref = [];

    frame.Omega = 0;
    frame.accel = 0;

    mirror(1).loc = 0.15;
    mirror(2).loc = -0.15;
    for k = 1:2
        mirror(k).dtheta_x = 0;
        mirror(k).dtheta_y = 0;
        mirror(k).cmask = true(size(sim.X));
        mirror(k).rmask = ones(size(sim.X));
    end

    outputs.R1 = [];
    outputs.R2 = [];
    outputs.loss1 = [];
    outputs.loss2 = [];
    outputs.Imax = [];
    outputs.zs = [];
    outputs.ts = [];

    toggles.track_centers = false;
    toggles.outputs_switch = true;
    toggles.videoplot_frequency = 'never';
    toggles.absorbing_mask = false;
    toggles.gain_switch = false;

    [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, struct());
    verifyEqual(testCase, laser.pos, mirror(1).loc);
    verifyEqual(testCase, sim.direction, -1);

    [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, struct());
    verifyEqual(testCase, laser.pos, mirror(2).loc);
    verifyEqual(testCase, sim.direction, 1);
    path_length = mirror(1).loc + mirror(1).loc - mirror(2).loc;
    verifyEqual(testCase, laser.time, path_length/consts.c, 'AbsTol', 1e-15);
    verifyTrue(testCase, all(diff(outputs.ts) > 0));
    verifyEqual(testCase, numel(outputs.ts), numel(outputs.zs));
end

function testLegacyWorkspaceGetsTimeState(testCase)
    consts.c = 3e8;
    sim.dz = -0.12;
    sim.t0 = 2e-6;
    laser.pos = 0;
    frame = struct();
    mirror = struct();
    outputs = struct();
    gain_medium = struct();
    toggles = struct();

    workspace_path = [tempname '.mat'];
    cleanup = onCleanup(@() delete(workspace_path));
    save(workspace_path, 'consts', 'sim', 'laser', 'frame', 'mirror', ...
        'outputs', 'gain_medium', 'toggles');

    [~, sim, laser, ~, ~, outputs] = Load_Workspace(workspace_path);
    verifyEqual(testCase, sim.direction, -1);
    verifyEqual(testCase, sim.dt, abs(sim.dz)/consts.c);
    verifyEqual(testCase, laser.time, sim.t0);
    verifyEmpty(testCase, outputs.ts);
end
