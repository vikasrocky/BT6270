function hhmodel_spikes()
    current_range = 0:0.01:0.6;
    npeaks_hist = zeros(size(current_range));

    for cur_iter = 1:length(current_range)
        gkmax = 0.36;
        vk = -77;
        gnamax = 1.20;
        vna = 50;
        gl = 0.003;
        vl = -54.387;
        cm = 0.01;

        dt = 0.01;
        niter = 50000;
        t = (1:niter) * dt;

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
            gna = gnamax * m^3 * h;
            gk = gkmax * n^4;
            gtot = gna + gk + gl;
            vinf = ((gna * vna + gk * vk + gl * vl) + current_range(cur_iter)) / gtot;
            tauv = cm / gtot;

            v = vinf + (v - vinf) * exp(-dt / tauv);

            alpham = 0.1 * (v + 40) / (1 - exp(-(v + 40) / 10));
            betam = 4 * exp(-0.0556 * (v + 65));

            alphan = 0.01 * (v + 55) / (1 - exp(-(v + 55) / 10));
            betan = 0.125 * exp(-(v + 65) / 80);

            alphah = 0.07 * exp(-0.05 * (v + 65));
            betah = 1 / (1 + exp(-0.1 * (v + 35)));

            taum = 1 / (alpham + betam);
            tauh = 1 / (alphah + betah);
            taun = 1 / (alphan + betan);

            minf = alpham * taum;
            hinf = alphah * tauh;
            ninf = alphan * taun;

            m = minf + (m - minf) * exp(-dt / taum);
            h = hinf + (h - hinf) * exp(-dt / tauh);
            n = ninf + (n - ninf) * exp(-dt / taun);

            vhist(iteration) = v;
            mhist(iteration) = m;
            hhist(iteration) = h;
            nhist(iteration) = n;
        end

        % Counting peaks 
        npeaks = find_peaks(vhist);
        npeaks_hist(cur_iter) = npeaks;
    end

    [I1, I2, I3] = get_currents(npeaks_hist);
    fprintf('\nThe cutoff currents are as follows:\n');
    fprintf('I1 = %.2f microA/mm^2\n', current_range(I1));
    fprintf('I2 = %.2f microA/mm^2\n', current_range(I2));
    fprintf('I3 = %.2f microA/mm^2\n', current_range(I3));

    figure;
    hold on;
    plot([current_range(I1), current_range(I1)], [0, max(npeaks_hist)], 'r');
    plot([current_range(I2), current_range(I2)], [0, max(npeaks_hist)], 'r');
    plot([current_range(I3), current_range(I3)], [0, max(npeaks_hist)], 'r');
    plot(current_range, npeaks_hist);
    title('Frequency of Spiking vs Input Current');
    xlabel('Input Current (microA/mm^2)');
    ylabel('Frequency of Spiking');
    hold off;
end

function npeaks = find_peaks(vhist)
    threshold = 10;
    npeaks = 0;
    for i = 2:length(vhist)-1
        if (vhist(i) >= threshold) && (vhist(i) > vhist(i+1)) && (vhist(i) > vhist(i-1))
            npeaks = npeaks + 1;
        end
    end
end

function [I1, I2, I3] = get_currents(npeaks_hist)
    for i = 2:length(npeaks_hist)-1
        if npeaks_hist(i) > 0 && npeaks_hist(i-1) == 0
            I1 = i;
        end
        if npeaks_hist(i+1) - npeaks_hist(i) > 4
            I2 = i;
        end
        if npeaks_hist(i+1) - npeaks_hist(i) < -2
            I3 = i;
        end
    end
end
