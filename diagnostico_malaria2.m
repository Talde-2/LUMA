function [diagnostico, tipo, densidad] = diagnostico_malaria2(rutaImagen)
    % Leer la imagen desde la ruta proporcionada
    I = imread(rutaImagen);

    % Paso 1: Eliminar los glóbulos rojos
    I_gray = rgb2gray(I);
    I_contrast = adapthisteq(I_gray);
    th = graythresh(I_contrast);
    BW = imbinarize(I_contrast, th);
    BW = imopen(~BW, strel('disk', 1));
    BW_clean = bwareaopen(BW, 10);

    % Paso 2: Detección de parásitos
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

        circ = 0;
        if peri > 0
            circ = 4 * pi * area / (peri^2);
        end

        aspect_ratio = bbox(3) / bbox(4);
        if area > 80 && area < 600 && ...
                ecc < 0.9 && sol > 0.85 && ...
                circ > 0.6 && aspect_ratio > 0.9 && aspect_ratio < 1.35
            mask_parasitos(L == k) = true;
        end
    end

    % Calcular número de parásitos
    [~, num_parasitos] = bwlabel(mask_parasitos);

    % Clasificar diagnóstico
    if num_parasitos == 0
        diagnostico = 'Negativo';
        tipo = 'N/A';
        densidad = 0;
    else
        diagnostico = 'Positivo';

        % Clasificar tipo de malaria
        if num_parasitos < 10
            tipo = 'Plasmodium vivax';
        else
            tipo = 'Plasmodium falciparum';
        end

        % Calcular densidad parasitaria (ejemplo simplificado)
        densidad = num_parasitos * 100; % valor de ejemplo
    end
end