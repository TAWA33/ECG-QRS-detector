
% Pathologie 1
load("ecg_Path1.mat");
len1 = length(ecg);
t1 = linspace(0, len1/Fs, len1);

figure(1);
subplot(2,1,1);
plot(t1, ecg, 'b');
title('Pathologie 1');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;
xlim([0, t1(end)]);

% Pathologie 2
load("ecg_Path2.mat");
len2 = length(ecg);
t2 = linspace(0, len2/Fs, len2);

[f_dom, fc_bpm] = detect_fibrillation(ecg, Fs)

subplot(2,1,2);
plot(t2, ecg, 'r');
title('Pathologie 2 ');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;
xlim([0, t2(end)]);

sgtitle('Comparaison des deux types de fibrillation', 'FontSize', 15, 'FontWeight', 'bold');


function [freq_dominante_Hz, fc_bpm] = detect_fibrillation(ecg, Fs)
    N = length(ecg);
    f = (0:N/2-1) * Fs / N;
    spectre = abs(fft(ecg));
    spectre = spectre(1:N/2);
    
    % Bande de frequence de la FVentriculaire (240-600 bpm = 4-10 Hz)
    idx_bande = find(f >= 4 & f <= 10);
    
    % Trouver le pic maximal dans la bande
    [~, idx_max] = max(spectre(idx_bande));
    freq_dominante_Hz = f(idx_bande(idx_max));
    fc_bpm = freq_dominante_Hz * 60;
end


function delta_k = autocov_RR(RR_intervals, k)
    N = length(RR_intervals);
    if k < N
        mean_RR = mean(RR_intervals);
        somme = 0;
        for n = 1:(N - k)
            somme = somme + (RR_intervals(n+k) - mean_RR) * (RR_intervals(n) - mean_RR);
        end
        delta_k = somme / (N - k - 1);
    end
end