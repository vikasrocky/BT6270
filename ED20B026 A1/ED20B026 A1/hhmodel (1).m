current_range = [0.03, 0.06, 0.3, 0.46];

for cur_iter = 1:length(current_range)
    gkmax = 0.36;
    vk = -77;
    gnamax = 1.20;
    vna = 50;
    gl = 0.003;
    vl = -54.387;
    cm = 0.01;

    dt = 0.01;
    niter = 10000;
    t = (0:niter-1)*dt;

    v = -64.9964;
    m = 0.0530;
    h = 0.5960;
    n = 0.3177;

    gnahist = zeros(1, niter);
    gkhist = zeros(1, niter);
    vhist = zeros(1, niter);
    mhist = zeros(1, niter);
    hhist = zeros(1, niter);
    nhist = zeros(1, niter);

    for iteration = 1:niter
        gna = gnamax*m^3*h;
        gk = gkmax*n^4;
        gtot = gna + gk + gl;
        vinf = ((gna*vna + gk*vk + gl*vl) + current_range(cur_iter)) / gtot;
        tauv = cm / gtot;

        v = vinf + (v - vinf)*exp(-dt/tauv);

        alpham = 0.1*(v+40) / (1 - exp(-(v+40)/10));
        betam = 4*exp(-0.0556*(v+65));

        alphan = 0.01*(v+55) / (1 - exp(-(v+55)/10));
        betan = 0.125*exp(-(v+65)/80);

        alphah = 0.07*exp(-0.05*(v+65));
        betah = 1 / (1 + exp(-0.1*(v+35)));

        taum = 1 / (alpham + betam);
        tauh = 1 / (alphah + betah);
        taun = 1 / (alphan + betan);

        minf = alpham*taum;
        hinf = alphah*tauh;
        ninf = alphan*taun;

        m = minf + (m - minf)*exp(-dt/taum);
        h = hinf + (h - hinf)*exp(-dt/tauh);
        n = ninf + (n - ninf)*exp(-dt/taun);

        vhist(iteration) = v;
        mhist(iteration) = m;
        hhist(iteration) = h;
        nhist(iteration) = n;
    end

    figure;
    set(gca, 'FontSize', 16);  % Set font size for axis labels and tick labels
    plot(t, vhist);
    string = sprintf('Voltage variation vs time; I=%.2f microA/mm^2', current_range(cur_iter));
    title(string, 'FontSize', 16);  % Set font size for title
    xlabel('Time (ms)', 'FontSize', 16);  % Set font size for x-axis label
    ylabel('Voltage (mV)', 'FontSize', 16);  % Set font size for y-axis label

    figure;
    set(gca, 'FontSize', 16);
    plot(t, mhist,'y-');
    hold on;
    plot(t, hhist,'g.');
    plot(t, nhist,'b-');
    legend({'m', 'h', 'n'}, 'FontSize', 16);  % Set font size for legend
    string = sprintf('Gating variables vs time; I=%.2f microA/mm^2', current_range(cur_iter));
    title(string, 'FontSize', 16);
    xlabel('Time (ms)', 'FontSize', 16);
    ylabel('Gating variable probabilities', 'FontSize', 16);

    figure;
    set(gca, 'FontSize', 16);
    gna = gnamax*(mhist.^3).*hhist;
    gk = gkmax*nhist.^4;
    plot(t, gna);
    hold on;
    plot(t, gk);
    legend({'gna', 'gk'}, 'FontSize', 16);
    string = sprintf('Conductance variation vs time; I=%.2f microA/mm^2', current_range(cur_iter));
    title(string, 'FontSize', 16);
    xlabel('Time (ms)', 'FontSize', 16);
    ylabel('Conductance', 'FontSize', 16);
end