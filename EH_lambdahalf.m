%% ========================================================================
%  Feldsimulation eines vertikalen Lambda/2-Dipols über realem Erdboden
% ========================================================================
%
%  Beschreibung:
%  Dieses Skript berechnet die elektrische und magnetische Feldstärke
%  (|E| und |H|) eines vertikal aufgestellten Halbwellendipols oberhalb
%  eines verlustbehafteten Erdbodens.
%
%  Physikalisches Modell:
%  - Der Dipol wird in N kurze Stromsegmente unterteilt (Dipol-Segmentierung),
%    jedes Segment wird als infinitesimaler Hertzscher Dipol behandelt
%    (Momentenmethode / Dipol-Array-Ansatz).
%  - Die Stromverteilung entlang des Dipols folgt der klassischen
%    Sinus-Näherung für einen dünnen linearen Strahler.
%  - Der Einfluss des realen (verlustbehafteten) Bodens wird über die
%    Spiegelquellenmethode (Bildtheorie) berücksichtigt: Jedes reale
%    Stromsegment erhält ein Spiegelsegment unterhalb der Erdoberfläche.
%    Die Reflexion am Boden wird über den Fresnel-Reflexionskoeffizienten
%    für senkrechte (vertikale) Polarisation gewichtet.
%  - Zusätzlich wird die durch den Boden verursachte Änderung der
%    Eingangsimpedanz (Strahlungswiderstand + Blindanteil) berechnet,
%    um daraus den tatsächlichen Speisestrom für eine vorgegebene
%    Sendeleistung Ptx zu bestimmen.
%
%  Ergebnis:
%  - Horizontale Feldkarten (|H|, |E|) in der Höhe der Dipolmitte
%  - Vertikale Feldkarten (|H|, |E|) in einer Ebene durch den Dipol
%  - Konturlinien für ausgewählte Feldstärkewerte (z. B. Grenzwerte
%    für Personenschutz) mit Markierung einer definierten Aufenthaltshöhe
%    (hier: z = 3 m über Grund, z. B. Kopfhöhe eines Menschen)
%
%  Zeitkonvention: exp(+j*w*t)  (technische/ingenieurmäßige Konvention)
%
%                                    DL8MCG Juli 2026  nicht verifiziert !
% ========================================================================

clear
clc

%% --------------------- Grundparameter der Anordnung --------------------
f       = 14e6;              % Betriebsfrequenz [Hz]  (14 MHz, Kurzwelle)
c       = 299792458;         % Lichtgeschwindigkeit im Vakuum [m/s]
lambda  = c/f;                % Wellenlänge [m]
k       = 2*pi/lambda;        % Kreiswellenzahl [1/m]
Ptx     = 100;                % Sendeleistung an der Antenne [W]

%% --------------------- Geometrie des Lambda/2-Dipols --------------------
% Der Dipol wird als dünner linearer Leiter entlang der z-Achse
% (vertikal) modelliert und in N diskrete Stromsegmente unterteilt.
L    = 10;                    % Gesamtlänge des Dipols [m] (~Lambda/2 bei 14 MHz)
z0   = 9;                     % Höhe der Dipolmitte über Grund [m]
N    = 201;                   % Anzahl der Stromsegmente (für Segmentierung)
zseg = linspace(z0-L/2, z0+L/2, N);  % z-Koordinaten der Segmentmittelpunkte [m]
dl   = L/(N-1);                % Länge eines einzelnen Stromsegments [m]

%% --------------------- Bodenparameter (elektrische Eigenschaften) --------
% Realer Erdboden als homogenes, verlustbehaftetes Medium beschrieben
% durch relative Permittivität und Leitfähigkeit (typische Werte für
% mittelmäßig leitenden Ackerboden).
epsr  = 13;                   % relative Permittivität des Bodens [-]
sigma = 5e-3;                 % Leitfähigkeit des Bodens [S/m]

%% --------------------- Stromverteilung (Referenzstrom I0 = 1 A) ---------
% Für einen dünnen linearen Dipol wird näherungsweise eine sinusförmige
% Stromverteilung entlang der Antenne angenommen (Standardnäherung der
% Antennentheorie für dünne Dipole).
I0_ref   = 1;                                              % Referenz-Speisestrom [A]
Iseg_ref = I0_ref*sin(k*(L/2-abs(zseg-z0)));                % Stromverteilung entlang des Dipols [A]

%% --------------------- Eingangsimpedanz inkl. Bodeneinfluss -------------
% Zself: Freiraum-Eingangsimpedanz eines Lambda/2-Dipols (Standardnäherung,
%        Strahlungswiderstand ca. 73 Ohm, Blindanteil ca. +42.5 Ohm für
%        einen idealisierten dünnen Halbwellendipol).
% dZ:    Zusätzliche Impedanzänderung durch den Einfluss des verlust-
%        behafteten Bodens, berechnet über die im Dipol influenzierte
%        Spannung aus dem reflektierten Feld der Spiegelquellen.
Zself = 73 + 1j*42.5;                                       % Freiraum-Näherung der Eingangsimpedanz [Ohm]
dZ    = deltaZ_ground(zseg, Iseg_ref, dl, f, epsr, sigma, I0_ref);  % Bodenbedingter Impedanzanteil [Ohm]
Zin   = Zself + dZ;                                          % Resultierende Eingangsimpedanz [Ohm]

%% --------------------- Tatsächlicher Speisestrom für Ptx ----------------
% Aus der gewünschten Sendeleistung Ptx und dem Realteil der
% Eingangsimpedanz (Strahlungswiderstand) wird der tatsächlich benötigte
% Speisestrom-Amplitudenwert I0 bestimmt (P = 0.5 * I0^2 * Re(Zin)).
I0   = sqrt(2*Ptx/real(Zin));                 % Tatsächlicher Speisestrom (Scheitelwert) [A]
Iseg = I0*sin(k*(L/2-abs(zseg-z0)));           % Skalierte Stromverteilung entlang des Dipols [A]

%% ========================================================================
%  Horizontale Feldkarte in Höhe der Dipolmitte (z = z0)
% ========================================================================
x = linspace(-5,5,101);   % x-Koordinaten des Auswertegitters [m]
y = linspace(-5,5,101);   % y-Koordinaten des Auswertegitters [m]

Hmap = zeros(length(y),length(x));   % Betrag des H-Feldes [A/m]
Emap = zeros(length(y),length(x));   % Betrag des E-Feldes [V/m]

% Punktweise Auswertung des Gesamtfeldes (Summe aller Dipolsegmente
% inkl. deren Spiegelsegmente) an jedem Gitterpunkt der Ebene.
for ix=1:length(x)
    for iy=1:length(y)
        r=sqrt(x(ix)^2+y(iy)^2);
        if r < 0.1
            % Auswertepunkte zu nahe an der Dipolachse werden
            % ausgeblendet (Nahfeld-Singularität der Punktquelle)
            Hmap(iy,ix)=NaN;
            continue
        end
        P=[x(ix),y(iy),z0];   % Beobachtungspunkt in der Horizontalebene
        [E,H]=dipoleFieldGround(P,zseg,Iseg,dl,f,epsr,sigma);
        Hmap(iy,ix)=norm(H);   % Betrag des komplexen H-Feldvektors
        Emap(iy,ix)=norm(E);   % Betrag des komplexen E-Feldvektors
    end
end

figure
contourf(x,y,Hmap,50)
colorbar
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
title('|H| in A/m   Horizontalebene z = 9m')

figure
contourf(x,y,Emap,50)
colorbar
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
title('|E| in V/m   Horizontalebene z = 9m')

%% ========================================================================
%  Vertikale Feldkarte (Schnittebene y = 0) über der Höhe z
% ========================================================================
% Diese Ebene zeigt den seitlichen Abstand x sowie die Höhe z über Grund
% und ist besonders relevant für die Beurteilung der Feldstärke in
% üblichen Aufenthaltshöhen von Personen (z. B. z = 3 m).
x = linspace(0.1,5,100);    % Seitlicher Abstand x von 0.1 m bis 5 m
z = linspace(0,18,150);     % Höhe über Grund z von 0 m bis 18 m

Hmap=zeros(length(z),length(x));
Emap=zeros(length(z),length(x));

for ix=1:length(x)
    for iz=1:length(z)
        P=[x(ix),0,z(iz)];   % Beobachtungspunkt in der Vertikalebene (y=0)
        [E,H]=dipoleFieldGround(P,zseg,Iseg,dl,f,epsr,sigma);
        Hmap(iz,ix)=norm(H);
        Emap(iz,ix)=norm(E);
    end
end

% --- Vollflächige Darstellung des H-Feldes ---
figure
contourf(x,z,Hmap,50)
colorbar
xlabel('Abstand seitlich x [m]')
ylabel('Höhe über Bodenmodell z [m]')
title('|H| in A/m  Vertikalebene')

% --- Konturlinie für einen definierten H-Feld-Grenzwert (0.1 A/m) ---
% Die zusätzliche horizontale Linie bei z = 3 m markiert eine typische
% Aufenthaltshöhe (z. B. Kopf-/Körperhöhe) zur Beurteilung, ob der
% Grenzwert dort über- oder unterschritten wird.
figure
contour(x,z,Hmap,[0.1 0.1],'LineWidth',2)
grid on
hold on
yline(3, 'r--', 'LineWidth', 1.5, 'Label', 'z = 3 m')
xlabel('Abstand seitlich x [m]')
ylabel('Höhe über Boden z [m]')
title('|H| in A/m  Vertikalebene')

% --- Vollflächige Darstellung des E-Feldes ---
figure
contourf(x,z,Emap,50)
colorbar
xlabel('Abstand seitlich x [m]')
ylabel('Höhe über Boden z [m]')
title('|E| in V/m  Vertikalebene')

% --- Konturlinie für einen definierten E-Feld-Grenzwert (28 V/m) ---
figure
contour(x,z,Emap,[28 28],'LineWidth',2)
grid on
hold on
yline(3, 'r--', 'LineWidth', 1.5, 'Label', 'z = 3 m')
xlabel('Abstand seitlich x [m]')
ylabel('Höhe über Boden z [m]')
title('|E| in V/m  Vertikalebene')


%% ========================================================================
%  Funktion: Bodenbedingte Änderung der Eingangsimpedanz
% ========================================================================
function dZ = deltaZ_ground(zseg, Iseg, dl, f, epsr, sigma, I0)
% deltaZ_ground  Berechnet die durch den verlustbehafteten Erdboden
% verursachte Änderung der Dipol-Eingangsimpedanz.
%
% Methode:
%   Für jedes reale Stromsegment m wird die vom (gespiegelten) Bild-
%   feld aller Segmente n induzierte E-Feldkomponente in z-Richtung
%   am Ort des Segments m aufsummiert. Über die Reziprozitäts-/
%   Induktionsmethode ergibt sich daraus die zusätzliche, durch den
%   Boden hervorgerufene Impedanz.
%
% Eingaben:
%   zseg  - z-Koordinaten der Dipolsegmente [m]
%   Iseg  - Stromverteilung entlang der Segmente [A]
%   dl    - Segmentlänge [m]
%   f     - Frequenz [Hz]
%   epsr  - relative Permittivität des Bodens [-]
%   sigma - Leitfähigkeit des Bodens [S/m]
%   I0    - Referenz-Speisestrom, auf den Iseg normiert ist [A]
%
% Ausgabe:
%   dZ    - Komplexe, bodenbedingte Impedanzänderung [Ohm]

    % Reflexionskoeffizient bei senkrechtem Einfall (theta = 0),
    % als Näherung für die Kopplung direkt unterhalb des Dipols.
    Gamma0 = reflectionVertical(0, f, epsr, sigma);

    N = length(zseg);
    Zsum = 0;
    for m = 1:N
        Pm = [0 0 zseg(m)];   % Ort des "Empfangs"-Segments m
        Ez_bild = 0;
        % Summation des von allen Spiegelsegmenten n am Ort m
        % erzeugten E-Feldes (z-Komponente)
        for n = 1:N
            Ps_img = [0 0 -zseg(n)];   % Position des Spiegelsegments n
            [En,~] = hertzDipoleSegment(Pm, Ps_img, Iseg(n), dl, f);
            Ez_bild = Ez_bild + En(3);
        end
        % Beitrag zum Impedanz-Summenterm (gewichtet mit Reflexions-
        % koeffizient und lokalem Segmentstrom, entspricht einer
        % induzierten EMK-Betrachtung)
        Zsum = Zsum + Gamma0*Ez_bild * Iseg(m) * dl;
    end
    dZ = -Zsum / I0^2;   % Normierung auf den Referenzstrom
end


%% ========================================================================
%  Funktion: Gesamtfeld eines Segment-Dipols über realem Boden
% ========================================================================
function [E,H] = dipoleFieldGround(P,zseg,Iseg,dl,f,epsr,sigma)
% dipoleFieldGround  Berechnet das komplexe E- und H-Feld am Punkt P,
% erzeugt durch den segmentierten Dipol UND dessen Spiegelbild im
% verlustbehafteten Erdboden (Bildtheorie mit Fresnel-Gewichtung).
%
% Für jedes reale Stromsegment wird zusätzlich ein Spiegelsegment
% unterhalb der Erdoberfläche berücksichtigt. Das vom Spiegelsegment
% ausgehende Feld wird mit dem winkelabhängigen Fresnel-Reflexions-
% koeffizienten (vertikale Polarisation) gewichtet, um die endliche
% Leitfähigkeit und Permittivität des Bodens korrekt abzubilden.
%
% Eingaben:
%   P     - Beobachtungspunkt [x y z] [m]
%   zseg  - z-Koordinaten der Dipolsegmente [m]
%   Iseg  - komplexe Stromverteilung entlang der Segmente [A]
%   dl    - Segmentlänge [m]
%   f     - Frequenz [Hz]
%   epsr  - relative Permittivität des Bodens [-]
%   sigma - Leitfähigkeit des Bodens [S/m]
%
% Ausgaben:
%   E - komplexer E-Feldvektor [Ex Ey Ez] [V/m]
%   H - komplexer H-Feldvektor [Hx Hy Hz] [A/m]

E = [0 0 0];
H = [0 0 0];
N = length(zseg);
for n = 1:N
    %% Beitrag des realen Stromsegments (direkte Welle)
    Ps = [0 0 zseg(n)];
    [En,Hn] = hertzDipoleSegment(P,Ps,Iseg(n),dl,f);
    E = E + En;
    H = H + Hn;

    %% Beitrag des Spiegelsegments (reflektierte Welle über Bildtheorie)
    Ps_img = [0 0 -zseg(n)];

    % Einfallswinkel der vom Spiegelsegment ausgehenden Welle relativ
    % zur Bodennormalen (z-Achse), benötigt für den Fresnel-Koeffizienten
    Rvec = P - Ps_img;
    r = norm(Rvec);
    costheta = abs(Rvec(3))/r;
    theta = acos(costheta);

    % Winkel- und materialabhängiger Fresnel-Reflexionskoeffizient
    % für vertikale (parallele) Polarisation
    Gamma = reflectionVertical(theta,f,epsr,sigma);
    [En,Hn] = hertzDipoleSegment(P,Ps_img,Iseg(n),dl,f);
    E = E + Gamma*En;   % Reflektiertes Feld wird mit Gamma gewichtet addiert
    H = H + Gamma*Hn;
end
end


%% ========================================================================
%  Funktion: Fresnel-Reflexionskoeffizient (vertikale Polarisation)
% ========================================================================
function Gamma = reflectionVertical(theta,f,epsr,sigma)
% reflectionVertical  Berechnet den komplexen Fresnel-Reflexionskoeffizienten
% für vertikal (parallel zur Einfallsebene) polarisierte Wellen an der
% Grenzfläche Luft/verlustbehafteter Boden.
%
% Der Boden wird über eine komplexe relative Permittivität beschrieben,
% die den ohmschen Verlustanteil (Leitfähigkeit sigma) mit einbezieht.
%
% Eingaben:
%   theta - Einfallswinkel gegenüber der Bodennormalen [rad]
%   f     - Frequenz [Hz]
%   epsr  - relative Permittivität des Bodens [-]
%   sigma - Leitfähigkeit des Bodens [S/m]
%
% Ausgabe:
%   Gamma - komplexer Reflexionskoeffizient [-]

eps0 = 8.854187817e-12;         % Elektrische Feldkonstante [F/m]
omega = 2*pi*f;                 % Kreisfrequenz [rad/s]

% Komplexe relative Permittivität des Bodens (Verlustanteil über
% Leitfähigkeit sigma gemäß epsc = epsr - j*sigma/(omega*eps0))
epsc = epsr - 1j*sigma/(omega*eps0);

s = sin(theta);
c = cos(theta);
root = sqrt(epsc - s.^2);

% Fresnel-Formel für parallele (vertikale) Polarisation
Gamma = (epsc*c - root) ./ ...
        (epsc*c + root);
end


%% ========================================================================
%  Funktion: Feld eines infinitesimalen Hertzschen Dipolsegments
% ========================================================================
function [E,H] = hertzDipoleSegment(P,Ps,I,dl,f)
% hertzDipoleSegment  Berechnet das exakte Nah-, Zwischen- und Fernfeld
% eines infinitesimal kurzen, vertikal orientierten Hertzschen Dipols.
%
% Es werden alle Feldterme (1/r, 1/r^2, 1/r^3) berücksichtigt, sodass
% die Funktion sowohl im Nahfeld als auch im Fernfeld korrekte
% Ergebnisse liefert.
%
% Zeitkonvention:
%       exp(+j*w*t)
%
% Dipolrichtung:
%       z-Richtung (vertikal)
%
% Eingaben:
%   P  = [x y z] Beobachtungspunkt       [m]
%   Ps = [x y z] Segmentposition         [m]
%   I  = komplexer Segmentstrom          [A]
%   dl = Segmentlänge                    [m]
%   f  = Frequenz                        [Hz]
%
% Ausgaben:
%   E = [Ex Ey Ez] komplexes E-Feld      [V/m]
%   H = [Hx Hy Hz] komplexes H-Feld      [A/m]

%% Naturkonstanten
mu0  = 4*pi*1e-7;              % Magnetische Feldkonstante [H/m]
eps0 = 8.854187817e-12;        % Elektrische Feldkonstante [F/m]
c    = 1/sqrt(mu0*eps0);       % Lichtgeschwindigkeit im Vakuum [m/s]
eta  = sqrt(mu0/eps0);         % Freiraum-Wellenwiderstand [Ohm] (~377 Ohm)

%% Wellenparameter
lambda = c/f;                  % Wellenlänge [m]
k = 2*pi/lambda;                % Kreiswellenzahl [1/m]

%% Abstand und Richtung vom Segment zum Beobachtungspunkt
Rvec = P - Ps;
r = norm(Rvec);
% Schutz gegen Division durch Null (Beobachtungspunkt im Segment)
if r < 1e-12
    error('Beobachtungspunkt liegt im Segment');
end

%% Kugelkoordinaten bezogen auf die Dipolachse (z-Richtung)
costheta = Rvec(3)/r;
sintheta = sqrt(1-costheta^2);

%% Feldkomponenten in Kugelkoordinaten (klassische Hertz-Dipol-Formeln)

% Azimutale Magnetfeldkomponente H_phi
% (Terme: Fernfeld 1/r  +  Nahfeld/Induktionsfeld 1/r^2)
Hphi = I*dl*sintheta/(4*pi) * ...
       (1j*k/r + 1/r^2) .* ...
       exp(-1j*k*r);

% Radiale elektrische Feldkomponente E_r
% (reines Nah-/Quasistatikfeld, fällt schneller ab als Fernfeld)
Er = eta * I*dl*costheta/(2*pi) * ...
     (1/r^2 - 1j/(k*r^3)) .* ...
     exp(-1j*k*r);

% Polare elektrische Feldkomponente E_theta
% (enthält Fernfeld-, Induktions- und Nahfeldanteil)
Etheta = 1j*eta*k*I*dl*sintheta/(4*pi) * ...
         (1/r - 1j/(k*r^2) - 1/(k^2*r^3)) .* ...
         exp(-1j*k*r);

%% Umrechnung von Kugel- in kartesische Koordinaten

% Radialer Einheitsvektor
er = Rvec/r;

% Einheitsvektor in theta-Richtung (tangential, in der Einfallsebene)
if sintheta < 1e-12
    % Auf der Dipolachse selbst (theta = 0) ist e_theta nicht definiert
    etheta = [0 0 0];
else
    etheta = [ ...
        costheta*Rvec(1)/(r*sintheta), ...
        costheta*Rvec(2)/(r*sintheta), ...
       -sintheta ];
end

% Einheitsvektor in phi-Richtung (azimutal, senkrecht zur Einfallsebene)
ephi = [-Rvec(2), Rvec(1), 0];
if norm(ephi) > 0
    ephi = ephi/norm(ephi);
end

%% Zusammensetzen der kartesischen Feldvektoren
Evec = Er*er + Etheta*etheta;
Hvec = Hphi*ephi;

%% Ausgabe
E = Evec;
H = Hvec;

end