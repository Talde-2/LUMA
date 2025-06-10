%% *SEGMENTACIÓN SIN GROUND TRUTH*

clear all 
close all

% Para optimizarlo y que vaya más rápido
set(0, 'DefaultFigureVisible', 'off'); % Desactivar visualización de figuras
warning('off', 'images:imhistc:inputHasNaNs'); % Desactivar warnings no críticos

files = dir(); % Coger todos los archivos
imageExtensions = {'.jpg', '.tiff'};
imageFiles = {};

for k = 1:length(files)
    [~, ~, ext] = fileparts(files(k).name);
    if ismember(lower(ext), imageExtensions)
        imageFiles{end+1} = files(k).name;
    end
end

columnas = {
    'ImageName','GlobulosBlancosDetectados', 'ParasitosDetectados', ... 
    'Leu_median_Area', 'Leu_median_BoundingBox','Leu_median_Eccentricity','Leu_median_Circularity','Leu_median_Diameter','Leu_median_Solidity','Leu_median_Perimeter','Leu_median_Intensity' ... 
    'Par_median_Area', 'Par_median_BoundingBox','Par_median_Eccentricity','Par_median_Circularity','Par_median_Diameter','Par_median_Solidity','Par_median_Perimeter','Par_median_Intensity'
};

Tabla_Global = table(...
    'Size', [0, numel(columnas)], ...
    'VariableTypes', [{'string'}, repmat({'double'}, 1, numel(columnas)-1)], ...
    'VariableNames', columnas ...
);
%% 
% Procesado de cada imagen:

for imgIdx = 1:length(imageFiles)
    I = imread(imageFiles{imgIdx});

    fprintf('\nProcesando imagen %d de %d: %s\n', imgIdx, length(imageFiles), imageFiles{imgIdx});

    % ============ 1. CONVERTIR A ESCALA DE GRISES ============ %
    [alto, ancho] = size(I);
    I_double = im2double(I);
    I_gris = rgb2gray(I);
    I_gris_double = im2double(I_gris);


    % ============ 2. APLICAR FILTROS ============ %
    I_Gauss = imgaussfilt(I_double);
    I_Gauss_double = im2double(I_Gauss);

    
    % ============ 3. OTSU ============ %
    A = mean(im2double(I_Gauss), 3); % Promedia los 3 canales RGB --> transforma a grises
    TO = graythresh(A);
    C = A > TO;

    % 3.1. Glóbulos Blancos
    BW_invert = ~C; % Invierte la imagen (blanco a negro y negro a blanco)
    BW_clean = bwareaopen(BW_invert, 50); % Quita ruido pequeño (con menos de 50 pixeles)
    [L, num] = bwlabel(BW_clean); % Etiqueta los objetos conectados = se identifican las zonas grandes
    
    stats = regionprops(L, 'Area', 'Centroid', 'BoundingBox');
    min_area = 100;
    max_area = 30000;
    mask_filtered = false(size(BW_clean));
    
    for k = 1:num
        if stats(k).Area > min_area && stats(k).Area < max_area
            mask_filtered(L == k) = true;
        end
    end % Nos quedamos con las formas que tienen un tamaño razoable = glóbulos blancos
    
    [L_filtered, num_filtered] = bwlabel(mask_filtered);
    stats_filtered = regionprops(L_filtered, 'Area', 'Centroid', 'BoundingBox');

    figure;
    imshow(I_Gauss); hold on;
    for k = 1:num_filtered
        rectangle('Position', stats_filtered(k).BoundingBox,'EdgeColor','g','LineWidth',2);
    end
    hold off;
    title(['Globulos blancos detectados:']);
    globulos_blancos_detectados = num2str(num_filtered)

    % 3.2. Parásitos
    [num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_Gauss)
    

    % ============ 4. GUARDAR CARACTERÍSTICAS ============ %
    stats_filtered_leu = regionprops(mask_filtered, I_gris_double, ...
    'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter', ...
    'Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');
    if ~isempty(stats_filtered_leu)
        Globulos_Blancos_Caracteristicas = struct2table(stats_filtered_leu);
    else
        Globulos_Blancos_Caracteristicas = table();
    end

    stats_filtered_par = regionprops(mask_parasitos, I_gris_double, ...
    'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter', ...
    'Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');
    if ~isempty(stats_filtered_par)
        Parasitos_Caracteristicas = struct2table(stats_filtered_par);
    else
        Parasitos_Caracteristicas = table();
    end

    nueva_fila = table(...
        string(imageFiles{imgIdx}), NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', columnas);
    nueva_fila.GlobulosBlancosDetectados = num_filtered;
    nueva_fila.ParasitosDetectados = num_parasitos;


    % 4.1. Glóbulos Blancos
    if ~isempty(stats_filtered_leu)
        nueva_fila.Leu_median_Area = median([stats_filtered_leu.Area]);
        nueva_fila.Leu_median_BoundingBox = mean(median([stats_filtered_leu.BoundingBox]));
        nueva_fila.Leu_median_Eccentricity = median([stats_filtered_leu.Eccentricity]);
        nueva_fila.Leu_median_Circularity = median([stats_filtered_leu.Circularity]);
        nueva_fila.Leu_median_Diameter = median([stats_filtered_leu.EquivDiameter]);
        nueva_fila.Leu_median_Solidity = median([stats_filtered_leu.Solidity]);
        nueva_fila.Leu_median_Perimeter = median([stats_filtered_leu.Perimeter]);
        nueva_fila.Leu_median_Intensity = median([stats_filtered_leu.MeanIntensity]);
    else
        nueva_fila{:, 4:11} = NaN;
    end
    
    % 4.2. Parásitos
    if ~isempty(stats_filtered_par)
        nueva_fila.Par_median_Area = median([stats_filtered_par.Area]);
        nueva_fila.Par_median_BoundingBox = mean(median([stats_filtered_par.BoundingBox]));
        nueva_fila.Par_median_Eccentricity = median([stats_filtered_par.Eccentricity]);
        nueva_fila.Par_median_Circularity = median([stats_filtered_par.Circularity]);
        nueva_fila.Par_median_Diameter = median([stats_filtered_par.EquivDiameter]);
        nueva_fila.Par_median_Solidity = median([stats_filtered_par.Solidity]);
        nueva_fila.Par_median_Perimeter = median([stats_filtered_par.Perimeter]);
        nueva_fila.Par_median_Intensity = median([stats_filtered_par.MeanIntensity]);
    else
        nueva_fila{:, 12:end} = NaN;
    end
    Tabla_Global = [Tabla_Global; nueva_fila];
end

% ============ ELIMINAR POSIBLES DUPLICADOS (POR SI ACASO) ============ %
[~, idx] = unique(Tabla_Global.ImageName, 'stable');
Tabla_Global = Tabla_Global(idx,:);

writetable(Tabla_Global, 'Falciparum_Data.csv');

fprintf('\nProcesamiento completado para todas las imágenes.\n');
display(Tabla_Global);
%%
% 3.2. Parásitos
function [num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_sin_globulos)
    % Convertir a escala de grises y mejorar contraste
    I_gray = rgb2gray(I_sin_globulos);
    I_contrast = adapthisteq(I_gray);
    th = graythresh(I_contrast);
    BW = imbinarize(I_contrast, th);
    BW = imopen(~BW, strel('disk', 1));
    BW_clean = bwareaopen(BW, 10);

    % Etiquetar y extraer características
    [L, num] = bwlabel(BW_clean);
    stats = regionprops(L, 'Area', 'Centroid', 'EquivDiameter', ...
                           'Eccentricity', 'Solidity', 'Perimeter', 'BoundingBox');

    mask_parasitos = false(size(BW_clean));

    for k = 1:num
        area = stats(k).Area;
        ecc = stats(k).Eccentricity;
        sol = stats(k).Solidity;
        peri = stats(k).Perimeter;
        bbox = stats(k).BoundingBox;

        % Circularidad: 1 = perfecto círculo
        circ = 0;
        if peri > 0
            circ = 4 * pi * area / (peri^2);
        end

        % Aspect ratio del bounding box
        aspect_ratio = bbox(3) / bbox(4); % ancho / alto

        % Filtro basado en características morfológicas
        if area > 80 && area < 600 && ...
           ecc < 0.9 && sol > 0.85 && ...
           circ > 0.6 && aspect_ratio > 0.9 && aspect_ratio < 1.35
            mask_parasitos(L == k) = true;
        end
    end

    [~, num_parasitos] = bwlabel(mask_parasitos);
    num_parasitos
    stats_filtered = regionprops(mask_parasitos, 'Centroid', 'EquivDiameter');
    centros = cat(1, stats_filtered.Centroid);
    radios = [stats_filtered.EquivDiameter] / 2;
end