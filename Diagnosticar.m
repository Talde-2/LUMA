function nueva_fila = Diagnosticar(rutaImagen)
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
A = mean(im2double(I_Gauss), 3); 
TO = graythresh(A);
C = A > TO;

% Glóbulos Blancos
BW_invert = ~C; 
BW_clean = bwareaopen(BW_invert, 50); 
[L, num] = bwlabel(BW_clean); 

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

% Parásitos
[num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_Gauss);

% Extraer características
stats_filtered_leu = regionprops(mask_filtered, I_gris_double, ...
    'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter', ...
    'Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');

stats_filtered_par = regionprops(mask_parasitos, I_gris_double, ...
    'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter', ...
    'Solidity', 'MeanIntensity', 'BoundingBox', 'Circularity');

% Inicializar fila de resultados
columnas = {
    'ImageName','GlobulosBlancosDetectados', 'ParasitosDetectados', ... 
    'Leu_median_Area', 'Leu_median_BoundingBox','Leu_median_Eccentricity','Leu_median_Circularity','Leu_median_Diameter','Leu_median_Solidity','Leu_median_Perimeter','Leu_median_Intensity' ... 
    'Par_median_Area', 'Par_median_BoundingBox','Par_median_Eccentricity','Par_median_Circularity','Par_median_Diameter','Par_median_Solidity','Par_median_Perimeter','Par_median_Intensity', ...
    'DiagnosticoBinario', 'TipoMalaria', 'DensidadParasitaria', 'CategoriaDensidad'
};

nueva_fila = table(...
    string(filename), num_filtered, num_parasitos, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    "", "", NaN, "", ...
    'VariableNames', columnas);

% Glóbulos blancos
if ~isempty(stats_filtered_leu)
   nueva_fila.Leu_median_Area = median([stats_filtered_leu.Area]);
   nueva_fila.Leu_median_BoundingBox = mean(median([stats_filtered_leu.BoundingBox]));
   nueva_fila.Leu_median_Eccentricity = median([stats_filtered_leu.Eccentricity]);
   nueva_fila.Leu_median_Circularity = median([stats_filtered_leu.Circularity]);
   nueva_fila.Leu_median_Diameter = median([stats_filtered_leu.EquivDiameter]);
   nueva_fila.Leu_median_Solidity = median([stats_filtered_leu.Solidity]);
   nueva_fila.Leu_median_Perimeter = median([stats_filtered_leu.Perimeter]);
   nueva_fila.Leu_median_Intensity = median([stats_filtered_leu.MeanIntensity]);
end

% Parásitos
if ~isempty(stats_filtered_par)
   nueva_fila.Par_median_Area = median([stats_filtered_par.Area]);
   nueva_fila.Par_median_BoundingBox = mean(median([stats_filtered_par.BoundingBox]));
   nueva_fila.Par_median_Eccentricity = median([stats_filtered_par.Eccentricity]);
   nueva_fila.Par_median_Circularity = median([stats_filtered_par.Circularity]);
   nueva_fila.Par_median_Diameter = median([stats_filtered_par.EquivDiameter]);
   nueva_fila.Par_median_Solidity = median([stats_filtered_par.Solidity]);
   nueva_fila.Par_median_Perimeter = median([stats_filtered_par.Perimeter]);
   nueva_fila.Par_median_Intensity = median([stats_filtered_par.MeanIntensity]);
end

% ============ PREDICCIÓN ============ %
modelo_PBL_bin = load('model_LogisticRegression.mat');
modelo_PBL_bin = modelo_PBL_bin.model_LogisticRegression;

modelo_PBL = load('model_RandomForest100.mat');
modelo_PBL = modelo_PBL.model_RandomForest100;

x_nueva_bin = table2array(nueva_fila(:, {'Leu_median_Intensity','GlobulosBlancosDetectados','Leu_median_Circularity','Par_median_Intensity'}));
[predictions_bin, ~] = predict(modelo_PBL_bin, x_nueva_bin);
nueva_fila.DiagnosticoBinario = string(predictions_bin);

if predictions_bin == "enfermo"
    x_nueva = table2array(nueva_fila(:, {'GlobulosBlancosDetectados','Leu_median_Circularity','Par_median_Intensity'}));
    [predictions, ~] = predict(modelo_PBL.Trained{1}, x_nueva); 
    clase_final = string(predictions);
    nueva_fila.TipoMalaria = clase_final;

    % Densidad parasitaria
    parasitos = nueva_fila.ParasitosDetectados;
    leucocitos = nueva_fila.GlobulosBlancosDetectados;

    if (leucocitos + parasitos == 0)
       densidad_parasitaria = 0;
    else
       factor = 400;
       densidad_parasitaria = (parasitos / (leucocitos + parasitos)*factor);
    end
    nueva_fila.DensidadParasitaria = densidad_parasitaria;

    % Categoría según especie
    if strcmp(clase_final, 'falciparum')
       if densidad_parasitaria < 800
          categoria_densidad = 'Baja';
       elseif densidad_parasitaria <= 4000
          categoria_densidad = 'Moderada';
       else
          categoria_densidad = 'Alta';
       end
    elseif strcmp(clase_final, 'vivax')
       if densidad_parasitaria < 800
          categoria_densidad = 'Baja';
       elseif densidad_parasitaria <= 2400
          categoria_densidad = 'Moderada';
       else
          categoria_densidad = 'Alta';
       end
    else
       categoria_densidad = 'Desconocida';
    end
    nueva_fila.CategoriaDensidad = categoria_densidad;
else
    nueva_fila.TipoMalaria = "sano";
    nueva_fila.DensidadParasitaria = 0;
    nueva_fila.CategoriaDensidad = "sano";
end

end