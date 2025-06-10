function [diagnostico, tipo, densidad, imagen_segmentada] = diagnostico_malaria(rutaImagen)
%% *APLICACIÓN*
% *1. SEGMENTACIÓN DE IMAGEN* 

% ============ 1. CARGAR LA IMAGEN ============ %
filename = rutaImagen; 
I = imread(filename); 

% ============ 2. CONVERTIR A ESCALA DE GRISES ============ %
I_double = im2double(I);
I_gris = rgb2gray(I);
I_gris_double = im2double(I_gris);

% ============ 3. APLICAR FILTROS ============ %
I_Gauss = imgaussfilt(I_double);
I_Gauss_double = im2double(I_Gauss);

% ============ 4. OTSU ============ %
A = mean(im2double(I_Gauss), 3); % Promedia los 3 canales RGB --> transforma a grises
TO = graythresh(A);
C = A > TO;

% 4.1. Glóbulos Blancos
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
end

[L_filtered, num_filtered] = bwlabel(mask_filtered);
stats_filtered = regionprops(L_filtered, 'Area', 'Centroid', 'BoundingBox');

% 4.2. Parásitos
%[num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_Gauss)

I_sin_globulos = I_Gauss;
I_gray = rgb2gray(I_sin_globulos);
I_contrast = adapthisteq(I_gray);

th = graythresh(I_contrast);
BW = imbinarize(I_contrast, th);

BW = imopen(~BW, strel('disk', 1));
BW_clean = bwareaopen(BW, 10);

[L_par, num_parasitos] = bwlabel(BW_clean);
stats_parasitos = regionprops(L_par, 'Area', 'Centroid', 'EquivDiameter', ...
                       'Eccentricity', 'Solidity', 'Perimeter', 'BoundingBox');

mask_parasitos = false(size(BW_clean));

for k = 1:num_parasitos
    area = stats_parasitos(k).Area;
    ecc = stats_parasitos(k).Eccentricity;
    sol = stats_parasitos(k).Solidity;
    peri = stats_parasitos(k).Perimeter;
    bbox = stats_parasitos(k).BoundingBox;

    circ = 0;
    if peri > 0
        circ = 4 * pi * area / (peri^2);
    end

    aspect_ratio = bbox(3) / bbox(4);

    if area > 80 && area < 600 && ...
       ecc < 0.9 && sol > 0.85 && ...
       circ > 0.6 && aspect_ratio > 0.9 && aspect_ratio < 1.35
        mask_parasitos(L_par == k) = true;
    end
end
[~, num_parasitos] = bwlabel(mask_parasitos);
stats_filtered_par = regionprops(mask_parasitos, 'Centroid', 'EquivDiameter');
centros = cat(1, stats_filtered_par.Centroid);
radios = [stats_filtered_par.EquivDiameter] / 2;

imagen_segmentada = I; % imagen original

for i = 1:size(centros, 1)
    centro = centros(i, :);
    radio = radios(i);
    imagen_segmentada = insertShape(imagen_segmentada, ...
        'Circle', [centro, radio], ...
        'Color', 'red', 'LineWidth', 2);
end

% ============ 5. GUARDAR CARACTERÍSTICAS ============ %
stats_filtered_leu = regionprops(mask_filtered, I_gris_double,'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter','Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');
if ~isempty(stats_filtered_leu)
   Globulos_Blancos_Caracteristicas = struct2table(stats_filtered_leu);
else
   Globulos_Blancos_Caracteristicas = table();
end

stats_filtered_par = regionprops(mask_parasitos, I_gris_double,'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter','Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');
if ~isempty(stats_filtered_par)
   Parasitos_Caracteristicas = struct2table(stats_filtered_par);
else
   Parasitos_Caracteristicas = table();
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

nueva_fila = table(...
        string(filename), NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', columnas);

nueva_fila.GlobulosBlancosDetectados = num_filtered;
nueva_fila.ParasitosDetectados = num_parasitos;

% 5.1. Glóbulos blancos
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

% 5.2. Parásitos
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


% 2. DIAGNÓSTICO CON MODELOS PREENTRENADOS

% Cargar modelos solo una vez (fuera de este bloque si es posible)
modelo_PBL_bin = load('model_LogisticRegression.mat');
modelo_PBL_bin = modelo_PBL_bin.model_LogisticRegression;

%modelo_PBL = load('model_RandomForest100.mat','model_RandomForest100');
modelo_PBL = load('model_RandomForest100.mat');
modelo_PBL = modelo_PBL.model_RandomForest100;
disp(class(modelo_PBL));

modelo_PBL = modelo_PBL.Trained{1};

% Modelo binario
x_nueva_bin = table2array(nueva_fila(:, {'Leu_median_Intensity','GlobulosBlancosDetectados','Leu_median_Circularity','Par_median_Intensity'}));
[predictions_bin, scores_bin] = predict(modelo_PBL_bin, x_nueva_bin);

% Random Forest
x_nueva = table2array(nueva_fila(:, {'GlobulosBlancosDetectados','Leu_median_Circularity','Par_median_Intensity'}));
[predictions, scores] = predict(modelo_PBL, x_nueva); 

% Guardar predicciones en la fila
nueva_fila.DiagnosticoBinario = string(predictions_bin);
fprintf('Predicción binaria: %s\n', nueva_fila.DiagnosticoBinario);

if predictions_bin == 'enfermo'
    
    clase_final = string(predictions);
    nueva_fila.TipoMalaria = clase_final;

    % Densidad parasitaria
    parasitos = nueva_fila.ParasitosDetectados;
    leucocitos = nueva_fila.GlobulosBlancosDetectados;

    if (leucocitos + parasitos == 0)
       densidad_parasitaria = 0;
    else
       factor = 400; % es el mas comun --> convierte la proporcion de parasitos en una densidad parasitaria
       densidad_parasitaria = (parasitos / (leucocitos + parasitos)*factor);
    end
    nueva_fila.DensidadParasitaria = densidad_parasitaria;

    if strcmp(clase_final, 'falciparum')
       if densidad_parasitaria < 800
          categoria_densidad = 'Baja';
       elseif densidad_parasitaria <= 4000
          categoria_densidad = 'Moderada';
       else
          categoria_densidad = 'Alta';
       end
    end
    if strcmp(clase_final, 'vivax')
       if densidad_parasitaria < 800
          categoria_densidad = 'Baja';
       elseif densidad_parasitaria <= 2400
          categoria_densidad = 'Moderada';
       else
          categoria_densidad = 'Alta';
       end
    end
    nueva_fila.CategoriaDensidad = categoria_densidad;
    
        
    fprintf('Tipo de malaria: %s\n', clase_final);
    fprintf('Densidad parasitaria: %.2f parásitos/µL\n', densidad_parasitaria);
    fprintf('Categoría de densidad: %s\n', categoria_densidad);
else
    clase_final = "sano";
    densidad_parasitaria = "sano";
    nueva_fila.TipoMalaria = "sano";
    nueva_fila.DensidadParasitaria = 0;
    nueva_fila.CategoriaDensidad = "sano";
    fprintf('Diagnóstico final: sano\n');
end
diagnostico = predictions_bin;
tipo = clase_final;
densidad=densidad_parasitaria;
%densidad = sprintf('%.2f parásitos/µL', densidad_parasitaria);
