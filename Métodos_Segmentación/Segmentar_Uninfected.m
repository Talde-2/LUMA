clear all 
close all

% Para optimizarlo y que vaya más rápido
set(0, 'DefaultFigureVisible', 'off'); % Desactivar visualización de figuras
warning('off', 'images:imhistc:inputHasNaNs'); % Desactivar warnings no críticos

files = dir(); % Coger todos los archivos

imageExtensions = {'.jpg', '.tiff'};
textExtension = '.txt';

imageFiles = {};
textFiles = {};

columnas = {
    'ImageName','GlobulosBlancosDetectados', 'ParasitosDetectados', ... 
    'Leu_median_Area', 'Leu_median_BoundingBox','Leu_median_Eccentricity','Leu_median_Circularity','Leu_median_Diameter','Leu_median_Solidity','Leu_median_Perimeter','Leu_median_Intensity' ... 
    'Par_median_Area', 'Par_median_BoundingBox','Par_median_Eccentricity','Par_median_Circularity','Par_median_Diameter','Par_median_Solidity','Par_median_Perimeter','Par_median_Intensity'
};

% Crear una tabla vacía con los nombres de las columnas
Tabla_Global = table('Size', [0, numel(columnas)], 'VariableTypes', repmat({'double'}, 1, numel(columnas)), 'VariableNames', columnas);
Tabla_Global.ImageName = string(Tabla_Global.ImageName); % Cambiar el tipo de ImageName a string

% Separar archivos de imagenes y de texto
for j = 1:length(files)
    [~, name, ext] = fileparts(files(j).name);
    if any(strcmpi(ext, imageExtensions))
        imageFiles{end+1} = files(j).name;
    elseif strcmpi(ext, textExtension)
        textFiles{end+1} = files(j).name;
    end
end
TP_total = 0;
FP_total = 0;
FN_total = 0;
GT_total = 0;
Procesado de cada imagen:
for imgIdx = 1:length(imageFiles)
    I = imread(imageFiles{imgIdx});

    par_validos = false(0);
    wbc_validos = false(0);
    
    fprintf('\nProcesando imagen %d de %d: %s\n', imgIdx, length(imageFiles), imageFiles{imgIdx});
    %disp(['Tamaño original de la imagen: ', num2str(size(I))]);

    % 1. Convertir a escala de grises
    [alto, ancho] = size(I);
    I_double = im2double(I);
    I_gris = rgb2gray(I);
    I_gris_double = im2double(I_gris);

    % 2. Aplicar filtros
    I_Gauss = imgaussfilt(I_double);
    I_Gauss_double = im2double(I_Gauss);
    figure;
    %subplot(1,2,1); imshow(I); title(['Original: ', imageFiles{imgIdx}]);
    %subplot(1,2,2); imshow(I_Gauss); title('Gauss');

    % Buscar el archivo de texto correspondiente
    [~, baseName, ~] = fileparts(imageFiles{imgIdx});
%     correspondingTextFile = [baseName, textExtension];
    baseNameClean = strrep(baseName, '_', ''); % para que si no tienen _ tambien los coja
    correspondingTextFile = '';
    for t = 1:length(textFiles)
        [~, textBaseName, ~] = fileparts(textFiles{t});
        textBaseNameClean = strrep(textBaseName, '_', '');
        if strcmpi(baseNameClean, textBaseNameClean)
            correspondingTextFile = textFiles{t};
            break;
        end
    end


     %if ismember(correspondingTextFile, textFiles)
     if ~isempty(correspondingTextFile)
        % Leer el archivo de texto
        textContent = fileread(correspondingTextFile);
        %fprintf('Archivo de texto encontrado: %s\n', correspondingTextFile);
        

        [alto, ancho] = size(I_Gauss);
        lines = splitlines(textContent); 
        % Tabla 1
        lines1 = lines(1);
        dims = str2double(split(lines1, ','));
        table_dims = array2table(dims');
        % Tabla 2 
        parasite_lines = lines(contains(lines, 'Parasite') | contains(lines, 'Parasitized'));
        if isempty (parasite_lines)
            parasite_data = cell(0,9);
            num_parasitos = 0;
        else
            parasite_data = cell(length(parasite_lines), 9);
            for parIdx = 1:length(parasite_lines)
                parts = strsplit(parasite_lines{parIdx}, ',');
                parasite_data(parIdx,:) = parts(1:9);
            end
        end     
        parasite_table = cell2table(parasite_data);
        parasite_table.Properties.VariableNames = {'Tipo', 'Area', 'BBoxX', 'BBoxY', 'BBoxW', 'CentroX', 'CentroY', 'PerimetroX', 'PerimetroY'};

        % Tabla 3
        wbc_lines = lines(contains(lines, 'White_Blood_Cell'));
        wbc_data = cell(length(wbc_lines), 7);
        for wbcIdx = 1:length(wbc_lines)
            parts = strsplit(wbc_lines{wbcIdx}, ',');
            wbc_data(wbcIdx,:) = parts(1:7);
        end
        wbc_table = cell2table(wbc_data);

        % Parasitos marcados
        Parasitos_Anotaciones_Figura = figure;
        imshow(I_Gauss);
        hold on;
        x_centers = str2double(parasite_data(:,6));
        y_centers = str2double(parasite_data(:,7));
        x_perim = str2double(parasite_data(:,8));
        y_perim = str2double(parasite_data(:,9));
        radios = sqrt((x_perim - x_centers).^2 + (y_perim - y_centers).^2);
        for i = 1:length(x_centers)
            % Dibujar círculo con centro (x_centers(i), y_centers(i)) y radio radios(i)
            viscircles([x_centers(i), y_centers(i)], radios(i), 'Color', 'blue', 'LineWidth', 1);
            
            % Opcional: marcar el centro con un punto
            %plot(x_centers(i), y_centers(i), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
            
            % Opcional: marcar el punto del perímetro
            %plot(x_perim(i), y_perim(i), 'bo', 'MarkerSize', 8, 'LineWidth', 1);
        end
        hold off;
        title('Parásitos marcados con círculos azules');
        % White Blood Cells marcados
        % Leucocitos_Anotaciones_Figura = figure;
        % imshow(I_Gauss);
        % hold on; 
        % x_coords = str2double(wbc_table{:,6});
        % y_coords = str2double(wbc_table{:,7});
        % plot(x_coords, y_coords, 'r.', 'MarkerSize', 10);
        % hold off;
        % title('Glóbulos blancos marcados con puntos rojos');

        [num_globulos_blancos, mask_globulos, I_sin_globulos, centroides_detectados] = detectarGlobulosBlancos(I);
        [num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_sin_globulos);
        % Convertir coordenadas del ground truth a número
        x_gt = str2double(parasite_data(:,6));
        y_gt = str2double(parasite_data(:,7));
        gt_centros = [x_gt, y_gt];
        
        % Tolerancia para considerar una detección como acierto
        tolerancia = 25;

        % Marcar coincidencias
        TP = 0;
        FN = 0;
        FP = 0;
        
        gt_usado = false(size(gt_centros, 1), 1);  % para evitar contar múltiples veces un GT
        
        for i = 1:size(centros, 1)
            distancias = sqrt((gt_centros(:,1) - centros(i,1)).^2 + (gt_centros(:,2) - centros(i,2)).^2);
            [min_dist, idx] = min(distancias);
        
            if (min_dist < tolerancia) & (~gt_usado(idx))
                TP = TP + 1;
                gt_usado(idx) = true;
            else
                FP = FP + 1;  % No coincide con ningún GT válido
            end
        end
        
        FN = sum(~gt_usado);  % Ground truths no detectados
        
        % Cálculo de métricas
        Precision = TP / (TP + FP);
        Recall = TP / (TP + FN);
        % 
        fprintf('\n--- Validación para %s ---\n', imageFiles{imgIdx});
        fprintf('True Positives: %d\n', TP);
        fprintf('False Positives: %d\n', FP);
        fprintf('False Negatives: %d\n', FN);
        fprintf('Precision: %.2f\n', Precision);
        fprintf('Recall: %.2f\n', Recall);
        else
        fprintf('No se encontró archivo de texto para: %s\n', imageFiles{imgIdx});
     end
end
% Tp=mean(TP);
% Fp=mean(FP);
% Fn=mean(FN);
% precision=mean(Precision);
% recall=mean(Recall);
% fprintf('True Positives: %d\n', Tp);
% fprintf('False Positives: %d\n', Fp);
% fprintf('False Negatives: %d\n', Fn);
% fprintf('Precision: %.2f\n', precision);
% fprintf('Recall: %.2f\n', recall);


function [num_globulos_blancos, mask_globulos, I_sin_globulos, centroides_detectados] = detectarGlobulosBlancos(I_Gauss)
    % Convertir a escala de grises promedio
    A = mean(im2double(I_Gauss), 3);
    TO = graythresh(A);
    C = A > TO;

    % Invertir y limpiar
    BW_invert = ~C;
    BW_clean = bwareaopen(BW_invert, 50);

    % Etiquetar y filtrar regiones por tamaño
    [L, num] = bwlabel(BW_clean);
    stats = regionprops(L, 'Area', 'Centroid', 'BoundingBox');
    areas = [stats.Area];
    min_area = prctile(areas, 10);
    max_area = prctile(areas, 90);

    mask_globulos = false(size(BW_clean));
    centroides_detectados = [];  % ← Aquí se guardarán los centróides

    for k = 1:num
        if stats(k).Area > min_area && stats(k).Area < max_area
            mask_globulos(L == k) = true;
            centroides_detectados = [centroides_detectados; stats(k).Centroid];  % ← Guardar centróide
        end
    end

    [~, num_globulos_blancos] = bwlabel(mask_globulos);

    % Dilatar y eliminar de la imagen
    se = strel('disk', 20);
    mascara_dilatada = imdilate(mask_globulos, se);
    I_sin_globulos = I_Gauss;

    for c = 1:3
        canal = I_Gauss(:,:,c);
        I_sin_globulos(:,:,c) = regionfill(canal, mascara_dilatada);
    end
end

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
    
    stats_filtered = regionprops(mask_parasitos, 'Centroid', 'EquivDiameter');
    centros = cat(1, stats_filtered.Centroid);
    radios = [stats_filtered.EquivDiameter] / 2;
end

