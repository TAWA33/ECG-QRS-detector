load("ecg_Path2.mat");
set(gca, 'Color', 'none');  
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'InvertHardCopy', 'off');
%ecg=-ecg;

duree_fenetre=4; %sec
N=Fs*duree_fenetre;
w=hamming(N);
d=round(N/2);
Nfft=1024;



%ecg = filtre_encoche(ecg(1, 1:end), 0.95, 12.5, Fs);

[Sx,f,t]=spectro(ecg, Fs, w, d, Nfft);

figure;
imagesc(t, f, log10(Sx)); % En base logarithmique
axis xy;
xlabel('Temps (s)');
ylabel('Frequence (Hz)');
title('Spectrogramme ECG');
colorbar;


function [Sx,f,t] = spectro(x, Fs, w, d, Nfft)
    % x : signal entre ; fen : fenetre ustilisee ; Fs : frequence
    % echantillonage ; d : decalage entre les fenetres ; Nfft : nombre de
    % points en frequence
    % Sx matrice du spectogramme ((Nfft/2) + 1) * M ; f : frequence ; t :
    % temps au debut de chaque fenetre
    
    N=length(w); % Taille de la fenetre
    L=length(x); % Taille du signal d'entree
    M=length((1:d:L-N)); % Nombre de fenetres glissantes
    y=zeros(N,M); % M Trames non fenetrees de taille N
    
    % Etape 1
    % Decomposition du vecteur x
    for i=0:M-1
        y(:,i+1)=x(i*d+1 : i*d+1 + N-1)';
    end
    
    % Etape 2
    % Fenetrage de tous les y(:,i)
    for i=1:M
        y(:,i) = (y(:,i) .* w);
    end
    % Etape 3
    % Calcul de la TFD
    Y=zeros(Nfft,M);
    for i=1:M
        Y(:,i)=fft(y(:,i), Nfft);
    end
    f=(1:floor(Nfft/2)).*(Fs/Nfft); % Axe frequenciel
    Sx=(1/Nfft)*(abs(Y).^2); % Calcul de la DSP
    Sx=Sx(1:floor(Nfft/2),:); % Selection de l'intervale de frequence [0;Fs/2]
    t=(1:d:L-N)./Fs; % Axe temporel
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