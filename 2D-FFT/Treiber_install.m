%% installIAAdapters.m
% Lädt die Support-Pakete für 'winvideo' und 'gentl' herunter und installiert sie.

function Treiber_install()
    % Liste der Adapter mit eindeutigen Dateinamen und aktuellen URLs
    adapters = { ...
        struct( ...
            'Name','winvideo', ...
            'Pkg','OS Generic Video Interface', ...
            'URL','https://ssd.mathworks.com/supportfiles/imageacquisition/winvideo/install.mlpkginstall' ...
        ), ...
        struct( ...
            'Name','gentl', ...
            'Pkg','GenICam Interface', ...
            'URL','https://ssd.mathworks.com/supportfiles/imageacquisition/gentl/install.mlpkginstall' ...
        ) ...
    };

    % Eingebaute Adapter abfragen
    hw     = imaqhwinfo;
    have   = hw.InstalledAdaptors;

    % Temp-Ordner anlegen
    tmp = fullfile(tempdir,'ia_adapters');
    if ~exist(tmp,'dir'), mkdir(tmp); end

    % Schleife über alle Adapter
    for k = 1:numel(adapters)
        ad = adapters{k};

        if ismember(ad.Name, have)
            fprintf('✔ Adapter "%s" ist bereits installiert.\n', ad.Name);
            continue
        end

        fprintf('– Installiere "%s" (%s)…\n', ad.Pkg, ad.Name);

        % Eindeutigen Namen für die Installationsdatei anlegen
        [~,~,ext] = fileparts(ad.URL);
        pkgFile   = fullfile(tmp, sprintf('%s_install%s', ad.Name, ext));

        % herunterladen, falls nicht vorhanden
        if ~exist(pkgFile,'file')
            try
                fprintf('  Download von:\n    %s\n', ad.URL);
                websave(pkgFile, ad.URL);
                fprintf('  Gespeichert als: %s\n', pkgFile);
            catch ME
                warning('MyScript:InstallError', '%s: %s', ME.identifier, ME.message);
            end

        else
            fprintf('  Nutze vorhandene Datei: %s\n', pkgFile);
        end

        % Installation versuchen (R2022b+ silent API)
        if exist('matlab.addons.supportpackage.installSupportPackage','file')
            try
                inst = matlab.addons.supportpackage.installSupportPackage(pkgFile);
                inst.waitForCompletion(300);  % Warte bis zu 5 Minuten
                fprintf('  Silent-Installation abgeschlossen.\n');
            catch ME
                
                supportPackageInstaller(pkgFile);
                promptContinue(ad);
            end
        else
            % Fallback: Support-Package Installer UI ohne Argument
            fprintf('  Starte Support-Package-Installer …\n');
            try
                % Versuch – neuer API nimmt evtl. pkgFile, ansonsten Fehler
                supportPackageInstaller(pkgFile);
            catch
                % Ältere MATLAB-Versionen: kein Arg erlaubt
                supportPackageInstaller;
            end

            % Hinweis: Im Installer "Install from Files" wählen und pkgFile auswählen
            promptContinue(ad);

        end
    end

    % Registry zurücksetzen und anzeigen
    imaqreset;
    final = imaqhwinfo;
    fprintf('\nInstallierte Adapter:\n');
    disp(final.InstalledAdaptors);
end

function promptContinue(ad)
    % Fordert Nutzer auf, den GUI-Dialog abzuschließen
    fprintf('\nBitte installiere "%s" im Support-Package-Installer und drücke ENTER.\n', ad.Pkg);
    pause;
end
