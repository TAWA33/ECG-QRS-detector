clear all; clc;

set(gca, 'Color', 'none');  
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'InvertHardCopy', 'off');

load("ecg_normal_1.mat");

ecg=-ecg;
fenetre_QRS=round(0.150*Fs);
ecg=ecg(1, 500:5000);
limite=length(ecg);
total_delay=23+ceil(fenetre_QRS/2);


%% Traitement du Signal pour la détection des impulsions R
x_enc = filtre_encoche(ecg(1, 1:end), 0.95, 12.5, Fs);
x_filt=pass_band_filter(ecg);
x_der=derivate(x_filt);
x_carr=x_der.^2;

x_integ=integration(x_carr, fenetre_QRS);

seuil_int=mean(x_integ)+0.5*sqrt(var(x_integ));

%Fonction avec seuil fixe
%imp_R = impulsion_R(ecg, x_integ, seuil_int, total_delay, fenetre_QRS);

%Fonction avec seuil adaptatif (SPKI / NPKI)
[imp_R, SPKI, NPKI, seuil_adaptatif] = impulsion_R_adaptatif(ecg, x_integ, total_delay, fenetre_QRS, Fs);

R_peaks = find(imp_R == 1);

%% Détection des impulsions Q et S à partir des R_peaks
[Q_peaks, S_peaks] = detect_Q_S(ecg, R_peaks, Fs);

%% Calcul de la fréquence cardiaque
calculate_heart_rate(R_peaks, Fs)

%% Calcul de la fonction d'autocovariance (Fibrillation auricuculaire Méthode 1)
RR_inter=diff(R_peaks)/Fs;
auto = zeros(1, length(RR_inter)-1);
auto0 = var(RR_inter); % Energie du signal

for i = 1:length(RR_inter)-1
    auto(i) = autocov_RR(RR_inter, i);
end

% Normalisation pour avoir le coefficient d'autocovariance
rho = auto / auto0;  % pour k=1..N-1
rho_complet = [1, rho];

%% Détection des impulsions P et T (Fibrillation auriculaire Méthode 2)

[P_peaks, T_peaks] = detect_P_T(ecg, R_peaks);
indices_P = (P_peaks > 0);
indices_T = (T_peaks > 0);

%% AFFICHAGE DES FIGURES
%==========================================

figure('Color', 'none');
plot(ecg);
hold on;
scatter(R_peaks, ecg(R_peaks), 'r^');
scatter(Q_peaks, ecg(Q_peaks), 'go');
scatter(S_peaks, ecg(S_peaks), 'bs');
scatter(P_peaks(indices_P), ecg(P_peaks(indices_P)), 'y');  % magenta pour P
scatter(T_peaks(indices_T), ecg(T_peaks(indices_T)), 'c');  % cyan pour T
xlim([50 1000])
legend('ECG', 'R', 'Q', 'S', 'P', 'T');





figure;
plot(0:length(rho), rho_complet, 'x-');
xlabel('Décalage k');
ylabel('\rho_k');
title('Autocorrélation des intervalles RR (normalisée)');
ylim([-1 1]);
grid on;
%{

figure,

subplot(2,1,1)
plot(ecg(1,1:end))
title('Signal ECG original')
ylabel('Amplitude')
grid on


subplot(2,1,2)
plot(ecg(1,1:end))
hold on
scatter(R_peaks, ecg(R_peaks), 'r^')
scatter(Q_peaks, ecg(Q_peaks), 'go')
scatter(S_peaks, ecg(S_peaks), 'bs')
title('Détection des ondes Q, R et S')
ylabel('Amplitude')
legend('ECG', 'Onde R', 'Onde Q', 'Onde S', 'Location', 'best')
grid on


figure('Position', [50 50 1200 800])


subplot(4,1,1)
plot(ecg(1,1:end))
title('Signal ECG original')
ylabel('Amplitude')
grid on

subplot(4,1,2)
semilogy(x_integ(1:end))
hold on
yline(seuil_int, 'r', 'LineWidth', 2)
title('Signal intégré avec SEUIL FIXE')
ylabel('Amplitude')
legend('x\_integ', 'Seuil fixe')
grid on

subplot(4,1,3)
semilogy(x_integ(1:end))
hold on
%plot(seuil_adaptatif(1:end), 'r')
title('Signal intégré avec SEUIL ADAPTATIF SPKI/NPKI - Résiste à la dérive')
ylabel('Amplitude')
legend('x\_integ', 'Seuil adaptatif')
grid on

subplot(4,1,4)
plot(ecg(1,1:end))
hold on
indices_R = find(imp_R == 1);
indices_R_lim = indices_R(indices_R <= length(ecg));
scatter(indices_R_lim, ecg(indices_R_lim), 'rx')
title('Détection des pics R avec seuil adaptatif (robuste à la dérive lente)')
xlabel('Échantillons')
ylabel('Amplitude')
grid on
%}

%===============================================================
%===============================================================PLOTS
%{

figure

limite = 1000;
Fe = 100;
r = 0.95;
f0 = 12.5;

signal_original = ecg(1, 1:limite);
signal_filtre = filtre_encoche(signal_original, r, f0, Fe);

N = length(signal_original);
frequences = (-N/2:N/2-1) * (Fe/N);

FFT_original = fftshift(fft(signal_original));
FFT_filtre = fftshift(fft(signal_filtre));

subplot(2,1,1)
plot(frequences, abs(FFT_original))
title('Spectre du signal ECG original')
xlabel('Fréquence (Hz)')
ylabel('|FFT|')
grid on


subplot(2,1,2)
plot(frequences, abs(FFT_filtre))
title('Spectre du signal ECG après filtre encoche (12.5 Hz)')
xlabel('Fréquence (Hz)')
ylabel('|FFT|')
grid on
%================================================================
%{
figure
subplot(5,1,1)
plot(ecg(1,1:limite))
hold on
scatter(indices_R_limites, ecg(indices_R_limites), 'rx', 'LineWidth', 1.5)
title("Signal ECG avec pics R")
hold off
%}

subplot(5,1,2)
plot(x_filt(1,1:limite))
title("Signal filtré par un passe-bande")

subplot(5,1,3)
plot(x_der(1,1:limite))
title("Signal dérivé")

subplot(5,1,4)
plot(x_carr(1,1:limite))
title("Signal dérivé au carré")

subplot(5,1,5), hold on
plot(x_integ(1,1:limite))
plot(1:0.7*10^3, ones(1,0.7*10^3)*seuil_int)
title("Signal intégré")


%saveas(gcf, "Localisation_RR_ECG.png")


%}

%==================================
%================================== FUNCTIONS


%% FILTRES
function x_filtered=pass_band_filter(x)
    Num_H_Low_Pass=zeros(1,13);
    Num_H_Low_Pass([1 7 13])=[1 -2 1];
    Den_H_Low_Pass=[1 -2 1];

    Num_H_High_Pass=zeros(1,33);
    Num_H_High_Pass([1 17 18 33])=[-1 32 -32 1];
    Den_H_High_Pass=[1 -1];

    x_inter=filter(Num_H_Low_Pass, Den_H_Low_Pass, x);
    x_filtered=filter(Num_H_High_Pass, Den_H_High_Pass, x_inter);

end

function x_deriv=derivate(x)
    Num_H_deriv=zeros(1,5);
    Num_H_deriv([1 2 4 5])=(1/8)*[-1 -2 2 1];

    x_deriv=filter(Num_H_deriv, 1, x);
end


function x_int = integration(x, N)
    x_int=zeros(1, length(x));
    for k=1:length(x)
        somme=0;
        compteur=0;
        for i=0:N-1
            if k-i >= 1
                somme=somme+x(k-i);
                compteur=compteur+1;
            end
        end
        x_int(k)=somme/compteur;
    end
end

function y = filtre_differentiateur(x)
    num = zeros(1, 7);
    num([1 7]) = [1 -1];
    denom = 1;
    y = filter(num, denom, x);
end

function y = filtre_passe_bas(x)
    num = [1, 0, 0, 0, 0, 0, 0, 0, -1];
    denom = [1 -1];
    y = filter(num, denom, x);
    
end

function x_enc = filtre_encoche(x, r, f, Fe)
    w0 = 2*pi*f/Fe;
    z1 = exp(1j*w0);
    z2 = exp(-1j*w0);
    p1 = r * exp(1j*w0);
    p2 = r * exp(-1j*w0);

    C = (1 - 2*r*cos(w0) + r^2) / (2 - 2*cos(w0));
    b = C * [1, -2*cos(w0), 1];
    a = [1, -2*r*cos(w0), r^2];

    x_enc = filter(b, a, x);
end


%% Détection des pics (Q, P, R, S, T)

% Fonction seuil Fixe
function imp_R = impulsion_R(x, x_integ, seuil_int, total_delay, fenetre_QRS)
    imp_R = zeros(1, length(x));
    
    indices_depassement = find(x_integ > seuil_int);
    
    zones_QRS = [];
    zone_actuelle = indices_depassement(1);
    debut_zone = zone_actuelle;
    
    for i = 2:length(indices_depassement)
        if indices_depassement(i) > indices_depassement(i-1) + fenetre_QRS/2

            zones_QRS = [zones_QRS; debut_zone, indices_depassement(i-1)];
            debut_zone = indices_depassement(i);
        end
    end
    zones_QRS = [zones_QRS; debut_zone, indices_depassement(end)];
    
    for k = 1:size(zones_QRS, 1)
        debut_zone = max(1, zones_QRS(k, 1) - total_delay);
        fin_zone = min(length(x), zones_QRS(k, 2) + total_delay);
        
        [~, idx_max] = max(abs(x(debut_zone:fin_zone)));
        idx_pic_R = debut_zone + idx_max - 1;
        
        imp_R(idx_pic_R) = 1;
    end
end

% Fonction seuil adaptatif
function [imp_R, SPKI, NPKI, thresholds_palier] = impulsion_R_adaptatif(x, x_integ, total_delay, fenetre_QRS, Fs)
    imp_R = zeros(1, length(x));
    
    decalage = 100;
    
    duree_init = min(round(1 * Fs), length(x_integ) - decalage);
    periode_init = decalage:(decalage + duree_init);
    
    [pics_init, ~] = findpeaks(x_integ(periode_init), 'MinPeakHeight', 0.2 * max(x_integ(periode_init)));
    
    if length(pics_init) >= 2
        SPKI = mean(pics_init);
    elseif isscalar(pics_init)
        SPKI = pics_init(1);
    else
        SPKI = max(x_integ(periode_init)) * 0.5;
    end
    
    NPKI = SPKI * 0.5;
    threshold = NPKI + 0.25 * (SPKI - NPKI);
    
    thresholds_palier = zeros(1, length(x_integ));
    thresholds_palier(1:decalage) = threshold;  % Remplir le debut
    
    [pics_valeurs, pics_indices] = findpeaks(x_integ, 'MinPeakDistance', round(0.2*Fs));
    
    dernier_idx = decalage;
    derniere_detection = -inf;
    periode_refractaire = round(fenetre_QRS * 0.5);
    
    for i = 1:length(pics_indices)
        idx_pic = pics_indices(i);
        valeur_pic = pics_valeurs(i);
        
        if idx_pic < decalage
            continue;
        end
        
        thresholds_palier(dernier_idx:idx_pic) = threshold;

        if (idx_pic - derniere_detection) < periode_refractaire
            continue;
        end
        
        if valeur_pic > threshold
            % C'est un QRS
            SPKI = 0.125 * valeur_pic + 0.875 * SPKI;

            debut_zone = max(1, idx_pic - total_delay);
            fin_zone = min(length(x), idx_pic + total_delay);
            [~, idx_max] = max(x(debut_zone:fin_zone));
            idx_vrai_R = debut_zone + idx_max -1;
            imp_R(idx_vrai_R) = 1;
            derniere_detection = idx_pic;
        else
            % C'est du bruit
            NPKI = 0.125 * valeur_pic + 0.875 * NPKI;
        end
        
        threshold = NPKI + 0.25 * (SPKI - NPKI);
        
        dernier_idx = idx_pic;
    end
    
    thresholds_palier(dernier_idx:end) = threshold;
end

function [Q_peaks, S_peaks] = detect_Q_S(ecg, R_peaks, Fs)
    for i = 1:length(R_peaks)
        debut_Q = max(1, R_peaks(i) - round(0.05*Fs));
        [~, idx_Q] = min(ecg(debut_Q:R_peaks(i)));
        Q_peaks(i) = debut_Q + idx_Q - 1;
        
        fin_S = min(length(ecg), R_peaks(i) + round(0.05*Fs));
        [~, idx_S] = min(ecg(R_peaks(i):fin_S));
        S_peaks(i) = R_peaks(i) + idx_S -1;
    end
end

function [P_peaks, T_peaks] = detect_P_T(ecg, R_peaks)
    x_diff = filtre_differentiateur(ecg);
    x_filtre = filtre_passe_bas(x_diff);
    
    % Ajustés selon votre observation
    delay_T = 44;
    delay_P = 55;
    
    nb_R = length(R_peaks);
    P_peaks = NaN(1, nb_R);
    T_peaks = NaN(1, nb_R);
    
    for i = 1:nb_R-1
        RR = R_peaks(i+1) - R_peaks(i);
        
        % Onde T
        debut_T = R_peaks(i);
        fin_T = min(R_peaks(i) + round(0.7 * RR), length(ecg));
        zone_T = x_filtre(debut_T:fin_T);
        passage_zero = find(zone_T(1:end-1) <= 0 & zone_T(2:end) > 0, 1, 'first');
        
        if ~isempty(passage_zero)
            idx_T = debut_T + passage_zero;
            idx_T = max(1, idx_T - delay_T+2);
            T_peaks(i) = idx_T;
        end
        
        % Onde P
        debut_P = fin_T + 1;
        fin_P = R_peaks(i+1);
        if debut_P <= fin_P && debut_P <= length(ecg)
            [~, idx_max] = max(x_filtre(debut_P:fin_P));
            idx_P = debut_P + idx_max - 1 + delay_P;
            P_peaks(i+1) = idx_P-4;
        end
    end
end

%% Fonction fréquence cardiaque

function heart_rate = calculate_heart_rate(R_peaks, Fs)
    RR_intervals = diff(R_peaks) / Fs;  % secondes
    mean_RR = mean(RR_intervals);
    heart_rate = 60 / mean_RR;  % bpm
end


%% Fonction d'autocovariance
function delta_k = autocov_RR(RR_intervals, k)
    N = length(RR_intervals);
    if k <= 0 || k >= N
        delta_k = NaN;
        return;
    end
    mean_RR = mean(RR_intervals);
    somme = 0;
    for n = 1:(N - k)
        somme = somme + (RR_intervals(n+k) - mean_RR) * (RR_intervals(n) - mean_RR);
    end
    delta_k = somme / (N - k - 1);
end







