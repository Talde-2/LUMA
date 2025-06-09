function Iout = otsu(Iin)
% === RUTAS === TF23_C108P69
imageFolder = 'C:\Users\Xabier\OneDrive\Escritorio\fotos\PF\Fotos\fotos';

% === ARCHIVOS DE IMAGEN ===
imageFiles = [dir(fullfile(imageFolder, '*.jpg')); dir(fullfile(imageFolder, '*.tiff')); dir(fullfile(imageFolder, '*.png'))];
TPs = [];
FPs = [];
Precisions = [];
Recalls = [];
FNR = [];
%% 
% *ANOTACIONES*

for i = 1:length(imageFiles)
    % Fitxategiaren izen osoa
    filename=fullfile(imageFolder, imageFiles(i).name);

    % Irudia kargatu
    I=Iin;
    I_Gauss=imgaussfilt(I);
    % figure;
    % imshow(I_Gauss); hold on;

    % Carpetas de las anotaciones
    [~, baseName, ~] = fileparts(imageFiles(i).name);
    ruta_GT = fullfile('C:\Users\Xabier\OneDrive\Escritorio\fotos\PF\Datos\datos', [baseName,'.txt']);
    
    % Datos de cada carpeta
    datos = readtable(ruta_GT, 'Delimiter', ',', 'HeaderLines', 1);

    % ===== PARÁSITOS =====
    solo_parasitos = datos(strcmp(datos.Var2, 'Parasite') | strcmp(datos.Var2, 'Parasitized'), :);
    % if ismember('Var8', datos.Properties.VariableNames) && ismember('Var9', datos.Properties.VariableNames)
    % x_parasitos = (solo_parasitos.Var6 + solo_parasitos.Var8) / 2;
    % y_parasitos = (solo_parasitos.Var7 + solo_parasitos.Var9) / 2;
    % plot(x_parasitos, y_parasitos, 'ro', 'MarkerSize', 5, 'LineWidth', 0.5); % círculos rojos
    % leyenda_parasitos = true;
    % else
    % disp(['No se encontraron coordenadas de parásitos en la imagen: ', imageFiles(i).name]);
    % leyenda_parasitos = false;
    % end
    
    %===== GLOBULOS BLANCOS =====
    solo_WBC = datos(strcmp(datos.Var2, 'White_Blood_Cell'), :);
    

    x_wbc = solo_WBC.Var6;
    y_wbc = solo_WBC.Var7;
    
    %plot(x_wbc, y_wbc, 'bo', 'MarkerSize', 5, 'LineWidth', 0.5);  % cruces azules
   

    % title(['Imagen: ', imageFiles(i).name], 'Interpreter', 'none');
    % if leyenda_parasitos
    % legend('Parásitos', 'WBC');
    % else
    % legend('WBC');
    % end
    % hold off;
   
    % ===== ESTRUCTURAS CELULARES =====
    [alto, ancho, ~] = size(I);  % Obtiene las dimensiones de la imagen original
    [Nucleos_Limpios, Citoplasma_Limpios, Nucleos, Citoplasma] = extraerEstructurasCelulares(I_Gauss, alto, ancho);
    % figure; --> funciona pero peta
    % subplot(2,2,1); imshow(Nucleos); title('Núcleos');
    % subplot(2,2,2); imshow(Nucleos_Limpios); title('Núcleos limpios');
    % subplot(2,2,3); imshow(Citoplasma); title('Citoplasma');
    % subplot(2,2,4); imshow(Citoplasma_Limpios); title('Citoplasma limpio');
    
    % Detectar glóbulos blancos y eliminar de la imagen
    [num_globulos_blancos, mask_globulos, I_sin_globulos, centroides_detectados] = detectarGlobulosBlancos(I);
    % %Parametros
    % tolerancia = 500;
    % 
    % % Extraer centróides ground truth (solo_WBC)
    % centroides_gt = [x_wbc, y_wbc];
    % 
    % % Inicializar confirmación
    % confirmados = [];
    % 
    % % Comparar cada GT con los detectados
    % for i = 1:size(centroides_gt, 1)
    %     gt = centroides_gt(i, :);
    %     distancias = vecnorm(centroides_detectados - gt, 2, 2);
    %     [min_dist, idx_min] = min(distancias);
    % 
    %     if min_dist < tolerancia
    %         confirmados = [confirmados; gt]; % Guardar centróide confirmado
    %     end
    % end
    % 
    % fprintf('Centróides confirmados: %d de %d\n', size(confirmados,1), size(centroides_gt,1));
    % 
    % % Visualización
    % figure; imshow(I_Gauss); hold on;
    % plot(centroides_detectados(:,1), centroides_detectados(:,2), 'bo', 'MarkerSize', 5);
    % plot(x_wbc, y_wbc, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5) % Ground truth
    % plot(confirmados(:,1), confirmados(:,2), 'go', 'MarkerSize', 8) % Confirmados
    % legend('Detecciones', 'GT', 'Confirmados');
    % title('Confirmación de centróides de glóbulos blancos');

    % Detectar parásitos Otsu
    [num_parasitos, mask_parasitos, centros, radios] = detectarParasitos(I_sin_globulos);
    % Detectar parásitos Canny
    %[num_parasitos, mask_parasitos, centros, radios] = detectarParasitosCanny(I_sin_globulos);
    % Detectar parásitos Gabor
    %[num_parasitos, mask_parasitos, centros, radios] = detectarParasitosGabor(I_sin_globulos);
    % Detectar parásitos k-means
    %[num_parasitos, mask_parasitos, centros, radios] = detectarParasitosKMeans(I_sin_globulos);

    % === Confirmar candidatos a parásitos con validación morfológica y de color ===
    % Confirmar candidatos comparando con GT (posición)
    % píxeles, ajustable
    %25 para otsu
    %50+ para Canny
    %80+ para K-means
    %105+ para Gabor
     tolerancia = 25;
     centros_confirmados = confirmarParasitosPorProximidad(centros, datos, tolerancia);
    figure; imshow(I_Gauss); hold on;
    plot(num_globulos_blancos, 'bo', 'MarkerSize', 5, 'LineWidth', 0.5);
    plot(centros_confirmados(:,1), centros_confirmados(:,2), 'r.', 'MarkerSize', 10);
    title(['Glóbulos blancos: ', num2str(num_globulos_blancos)-1, ...
           ' | Parásitos confirmados (GT): ', num2str(size(centros_confirmados,1))]);
    hold off;
    % total_wbc_gt = 0;
    % num_wbc_gt = height(solo_WBC);
    % total_wbc_gt = total_wbc_gt + num_wbc_gt;
    % fprintf('Imagen: %s | Glóbulos blancos (GT): %d\n', imageFiles(i).name, num_wbc_gt);
    % total_par = 0;
    % num_par = height(solo_parasitos);
    % total_par = total_par + num_par;
    % fprintf('Imagen: %s | Parasitos (GT): %d\n', imageFiles(i).name, num_par);
    % Mostrar resultados
    % figure; imshow(I); hold on;
    % plot(centros(:,1), centros(:,2), 'r.', 'MarkerSize', 10);
    % title(['Glóbulos blancos: ', num2str(num_globulos_blancos), ' |  Parásitos confirmados (GT): ', num2str(size(centros_confirmados,1))]);
    % hold off;
    solo_parasitos = datos(strcmp(datos.Var2, 'Parasite') | strcmp(datos.Var2, 'Parasitized'), :);
    x_gt = (solo_parasitos.Var6 + solo_parasitos.Var8) / 2;
    y_gt = (solo_parasitos.Var7 + solo_parasitos.Var9) / 2;
    centros_gt = [x_gt, y_gt];
    [TP_rate, FP_rate, precision, recall, TP, FP, FN, TPmean, FNRmean, FPmean, Precisionmean, Recallmean, FN_centros] = evaluarDeteccion(centros_confirmados, centros_gt, tolerancia);

    % Mostrar resultados por imagen
    % fprintf('TP Rate: %.2f%%\n', TPmean);
    % fprintf('FP Rate: %.2f%%\n', FPmean);
    % fprintf('Precisión: %.2f%%\n', Precisionmean);
    % fprintf('Recall: %.2f%%\n', Recallmean);
    % 
    % Acumular métricas
    
    TPs(end+1) = TPmean;
    FPs(end+1) = FPmean;
    Precisions(end+1) = Precisionmean;
    Recalls(end+1) = Recallmean;
    FNR(end+1) = FNRmean;
end
    % Calcular métricas promedio
    TPmean_total = mean(TPs);
    FPmean_total = mean(FPs);
    Precision_total = mean(Precisions);
    Recall_total = mean(Recalls);
    FNR_total = mean(FNR);
    
    
    fprintf('\n=== MÉTRICAS PROMEDIO GLOBALES ===\n');
    fprintf('TPmean: %.2f%%\n', TPmean_total);
    fprintf('FPmean: %.2f%%\n', FPmean_total);
    fprintf('Precisión media: %.2f%%\n', Precision_total);
    fprintf('Recall media: %.2f%%\n', Recall_total);
    fprintf('Falsos negativos media: %.2f%%\n', FNR_total)

%% 
% *EXTRAER ESTRUCTURAS CELULARES*

function [Nucleos_Limpios, Citoplasma_Limpios, Nucleos, Citoplasma] = extraerEstructurasCelulares(I_Gauss, alto, ancho)
    % Extrae núcleos y citoplasma de una imagen teñida con Giemsa

    % Separación de canales RGB
    R = I_Gauss(:,:,1);  
    G = I_Gauss(:,:,2);  
    B = I_Gauss(:,:,3);  

    % ===== NÚCLEOS =====
    Nucleos = 2*B - R - G;  % Enfatiza azul
    Nucleos = imadjust(Nucleos); 
    Nucleos_Binario = imbinarize(Nucleos, 'adaptive','Sensitivity', 0.3); 
    Nucleos_Limpios = imopen(Nucleos_Binario, strel('disk', 2));
    Nucleos_Limpios = imresize(Nucleos_Limpios, [alto, ancho]);

    % ===== CITOPLASMA =====
    Citoplasma = (R + 0.7*G) - 0.8*B;
    Citoplasma = imadjust(Citoplasma); 
    Citoplasma_Binario = imbinarize(Citoplasma, 0.08); 
    Citoplasma_Limpios = imclose(Citoplasma_Binario, strel('disk', 2));
    Citoplasma_Limpios = imresize(Citoplasma_Limpios, [alto, ancho]);
end
%% 
% *DETECTAR GLOBULOS BLANCOS*

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
%% 
% *1-METODO: OTSU*

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
%% 
% *2-METODO: CANNY*

function [num_parasitos, mask_parasitos, centros, radios] = detectarParasitosCanny(I_sin_globulos)
    % Detección de parásitos basada en bordes (Canny)

    % Convertir a escala de grises
    I_gray = rgb2gray(I_sin_globulos);
    
    % Filtrar con mediana para reducir ruido

    I_gauss = imgaussfilt(I_gray);
    I_filt = medfilt2(I_gauss,[3 3]);

    % Detectar bordes usando Canny
    BW_edges = edge(I_filt, 'Canny');

    % Cerrar bordes (morfología)
    BW_closed = imclose(BW_edges, strel('disk', 3));

 
    % Eliminar objetos pequeños
    BW_clean = bwareaopen(BW_closed, 300);

    % Etiquetar y filtrar regiones por tamaño
    [L, num] = bwlabel(BW_clean);
    stats = regionprops(L, 'Area', 'Centroid', 'EquivDiameter');
   areas = [stats.Area];
    min_area = prctile(areas, 10);
    max_area = prctile(areas, 90);


    mask_parasitos = false(size(BW_clean));
    centros = [];
    radios = [];

    for k = 1:num
        if stats(k).Area > min_area && stats(k).Area < max_area
            mask_parasitos(L == k) = true;
            centros = [centros; stats(k).Centroid];
            radios = [radios, stats(k).EquivDiameter / 2];
        end
    end

    num_parasitos = size(centros, 1);
end
%% 
% *3-METODO: K-MEANS*

function [num_parasitos, mask_parasitos, centros, radios] = detectarParasitosKMeans(I_sin_globulos)

    % === 1. Segmentación con K-means ===
    rng(1);  % Asegura resultados reproducibles
    k = 3;   % Número de clusters
    [cluster_idx, ~] = imsegkmeans(I_sin_globulos, k);
    
    % === 2. Crear máscara para cada cluster ===
    masks = false([size(cluster_idx), k]);
    for i = 1:k
        masks(:,:,i) = cluster_idx == i;
        
    end
    
    % === 3. Analizar cada cluster para ver cuál parece contener los parásitos ===
    % Criterios: áreas pequeñas, formas redondeadas
    min_area = 50;
    max_area = 800;  % puedes ajustarlo según el tamaño del parásito

    best_mask = false(size(cluster_idx));
   
    centros = [];
    radios = [];

    for i = 1:k
        % Limpiar máscara
        mask = bwareaopen(masks(:,:,i), min_area);
        
        % Etiquetar componentes conectados
        [L, num] = bwlabel(mask);
        stats = regionprops(L, 'Area', 'Centroid', 'EquivDiameter');
        
        for j = 1:num
            area = stats(j).Area;
            if area >= min_area && area <= max_area
                best_mask = best_mask | (L == j);
                centros = [centros; stats(j).Centroid];
                radios = [radios; stats(j).EquivDiameter / 2];
            end
        end
    end

    % === 4. Salidas ===
    mask_parasitos = best_mask;
    num_parasitos = size(centros, 1);

end

%% 
% *4-METODO: GABOR*

function [num_parasitos, mask_parasitos, centros, radios] = detectarParasitosGabor(I_sin_globulos)
    % Convertir a escala de grises
    Igray = rgb2gray(I_sin_globulos);
    Igray = im2double(Igray);

    % Crear banco de filtros Gabor
    wavelength = [2 4 8];     % Longitudes de onda
    orientation = 0:45:135;   % Orientaciones
    gaborBank = gabor(wavelength, orientation);

    % Aplicar filtros
    gaborMag = imgaborfilt(Igray, gaborBank);

    % Combinar respuestas: media de magnitudes
    gaborSum = mean(gaborMag, 3);

    % Normalizar
    gaborSum = mat2gray(gaborSum);

    % Segmentar con Otsu sobre la textura
    level = graythresh(gaborSum);
    mask = imbinarize(gaborSum, level);

    % Limpieza morfológica
    mask_clean = imopen(mask, strel('disk', 1));
    mask_clean = imclose(mask_clean, strel('disk', 2));
    mask_parasitos = bwareaopen(mask_clean, 30);
    mask_roi = crearMascaraCircular(size(mask_parasitos));  % Función auxiliar
    mask_parasitos = mask_parasitos & mask_roi;
    mask_parasitos = ~mask_parasitos;
    
    % Extraer características
    stats = regionprops(mask_parasitos, 'Centroid', 'EquivDiameter');

    centros = reshape([stats.Centroid], 2, []).';
    radios = [stats.EquivDiameter] / 2;
    num_parasitos = size(centros, 1);
end
function mask_circular = crearMascaraCircular(imSize)
    % imSize: tamaño de la imagen [filas, columnas]
    rows = imSize(1);
    cols = imSize(2);
    
    % Centro y radio
    [X, Y] = meshgrid(1:cols, 1:rows);
    centerX = round(cols / 2);
    centerY = round(rows / 2);
    radius = min(centerX, centerY) * 0.95;  % 95% del tamaño de la imagen

    % Crear máscara circular
    mask_circular = (X - centerX).^2 + (Y - centerY).^2 <= radius^2;
end

%% 
% *CONFIRMAR CANDIDATOS*

function centros_confirmados = confirmarParasitosPorProximidad(centros_detectados, datos_gt, tolerancia)
    % Filtrar solo las anotaciones de tipo "Parasite"
    solo_parasitos = datos_gt(strcmp(datos_gt.Var2, 'Parasite') | strcmp(datos_gt.Var2, 'Parasitized'), :);

    % Verificar que las columnas necesarias existan
    if ~ismember('Var6', datos_gt.Properties.VariableNames) || ...
       ~ismember('Var8', datos_gt.Properties.VariableNames) || ...
       ~ismember('Var7', datos_gt.Properties.VariableNames) || ...
       ~ismember('Var9', datos_gt.Properties.VariableNames)
        warning('Las columnas necesarias no existen en el GT');
        centros_confirmados = [];
        return;
    end

    % Calcular los centros reales del GT: promedio entre esquina superior izquierda y esquina inferior derecha
    x_gt = (solo_parasitos.Var6 + solo_parasitos.Var8) / 2;
    y_gt = (solo_parasitos.Var7 + solo_parasitos.Var9) / 2;
    centros_gt = [x_gt, y_gt];

    % Inicializar lista de centros confirmados
    centros_confirmados = [];

    for i = 1:size(centros_detectados, 1)
        centro = centros_detectados(i, :);

        % Calcular distancia euclidiana a cada centro del GT
        distancias = sqrt(sum((centros_gt - centro).^2, 2));

        % Confirmar si hay alguna anotación dentro del radio de tolerancia
        if any(distancias < tolerancia)
            centros_confirmados = [centros_confirmados; centro];
        end
    end
end
function [TP_rate, FP_rate, FNR, precision, recall, ...
          TP, FP, FN, ...
          TPmean, FPmean, FNRmean, Precisionmean, Recallmean, ...
          FN_centros] = evaluarDeteccion(centros_detectados, centros_gt, tolerancia)

    % Inicializar métricas
    TP = 0;
    FP = 0;
    FN = 0;
    FN_centros = [];

    % Si no hay detecciones ni GT, devolver ceros
    if isempty(centros_detectados) && isempty(centros_gt)
        TP_rate = 0;
        FP_rate = 0;
        FNR = 0;
        precision = 0;
        recall = 0;
        TPmean = 0;
        FPmean = 0;
        FNRmean = 0;
        Precisionmean = 0;
        Recallmean = 0;
        return;
    end

    % Comparar cada centro detectado con los GT
    centros_usados = false(size(centros_gt,1), 1);

    for i = 1:size(centros_detectados,1)
        distancias = vecnorm(centros_gt - centros_detectados(i,:), 2, 2);
        [min_dist, idx_min] = min(distancias);
        if min_dist < tolerancia && ~centros_usados(idx_min)
            TP = TP + 1;
            centros_usados(idx_min) = true;
        else
            FP = FP + 1;
        end
    end

    % Falsos negativos: centros_gt no detectados
    for i = 1:size(centros_gt,1)
        if ~centros_usados(i)
            FN = FN + 1;
            FN_centros = [FN_centros; centros_gt(i,:)]; % guardar coordenadas FN
        end
    end

    % Métricas base (evitar división por cero)
    precision = TP / max(TP + FP, 1);
    recall = TP / max(TP + FN, 1);
    TP_rate = TP / max(size(centros_gt,1), 1);
    FP_rate = FP / max(size(centros_detectados,1), 1);
    FNR = FN / max(size(centros_gt,1), 1);  % False Negative Rate

    % Métricas en porcentaje
    TPmean = TP_rate * 100;
    FPmean = FP_rate * 100;
    FNRmean = FNR * 100;
    Precisionmean = precision * 100;
    Recallmean = recall * 100;
end

Iout= I_Gauss;
end