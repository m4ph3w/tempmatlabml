clear;
clc;
close all;

%% Generate Rotated Rectangle Image

img_height = 100;
img_width = 150;
binaryImage = zeros(img_height, img_width);

% Rectangle dimensions
rect_height = 2;
rect_width = 10;

% Starting point (centered)
start_row = floor((img_height - rect_height) / 2) + 1;
start_col = floor((img_width - rect_width) / 2) + 1;

% Place the rectangle in the image
binaryImage(start_row:start_row+rect_height-1, start_col:start_col+rect_width-1) = 1;

% Rotate the image by 30 degrees
rotatedImage = imrotate(binaryImage, 30, 'crop');

% Display the rotated image
figure, imshow(rotatedImage);
title('Rotated Rectangle');

IMedian = logical(rotatedImage);
[m, n] = size(IMedian);

%% Finding Centroid

% Compute the area (M_00)
M_00 = sum(IMedian, "all");

% Expected area
expected_area = rect_height * rect_width;

% Compute the first moments M_10 and M_01
[i_indices, j_indices] = find(IMedian);

M_10 = sum(j_indices); % Sum over columns (x-coordinate)
M_01 = sum(i_indices); % Sum over rows (y-coordinate)

% Compute centroids
centroid_x = M_10 / M_00; % x-coordinate (column index)
centroid_y = M_01 / M_00; % y-coordinate (row index)

% Compute the second moments M_20, M_02, and M_11
M_20 = sum(j_indices.^2);
M_02 = sum(i_indices.^2);
M_11 = sum(i_indices .* j_indices);

% Compute central moments mu_20, mu_02, and mu_11
mu_20 = M_20 - centroid_x * M_10;
mu_02 = M_02 - centroid_y * M_01;
mu_11 = M_11 - centroid_x * M_01;

% Form the covariance matrix J
J = [mu_20, mu_11; mu_11, mu_02];

% Compute eigenvalues and eigenvectors
[V,D] = eig(J);

% Compute the lengths of the major and minor axes
a = 2 * sqrt(max(diag(D)) / M_00);
b = 2 * sqrt(min(diag(D)) / M_00);

% Compute the orientation angle
[~, idx] = max(diag(D)); % Find the index of the maximum eigenvalue
V_major = V(:, idx);
theta = atan2(V_major(2), V_major(1)); % Theta in radians

% Convert theta to degrees
orientation_angle = rad2deg(theta);
orientation_angle = mod(-orientation_angle, 180); % Adjust for image coordinates

% Now overlay the centroid on the original image
figure, imshow(rotatedImage);
hold on;
plot(centroid_x, centroid_y, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
text(centroid_x + 5, centroid_y, sprintf('(%.2f, %.2f)', centroid_x, centroid_y), 'Color', 'red');
hold off;
title('Rotated Rectangle with Centroid Marked');

%% Edge Detection and Circularity Calculation

% Use bwperim to get the perimeter pixels
perimeterImage = bwperim(IMedian);
figure, imshow(perimeterImage);
title('Perimeter Image');

% Calculate perimeter by counting the number of pixels in the perimeter
perimeter_pixels = sum(perimeterImage(:));

% Calculate the chain code and perimeter
chain_code = eight_neighbor_chain_code(perimeterImage);
perimeter_chain = calculate_perimeter_from_chain_code(chain_code);

% Expected perimeter
expected_perimeter = 2 * (rect_height + rect_width);

% Compute circularity
CIRCULARITY = (4 * pi * M_00) / (perimeter_chain^2);

% Display calculated values
fprintf('Calculated area (M_00): %f\n', M_00);
fprintf('Expected area: %f\n', expected_area);
fprintf('Calculated perimeter (from chain code): %f\n', perimeter_chain);
fprintf('Calculated perimeter (pixel count): %d\n', perimeter_pixels);
fprintf('Expected perimeter: %f\n', expected_perimeter);
fprintf('Circularity: %f\n', CIRCULARITY);
fprintf('Orientation angle (degrees): %f\n', orientation_angle);

%% Function Definitions

function chain_code = eight_neighbor_chain_code(binary_image)
    % Ensure the image is binary
    binary_image = logical(binary_image);

    % Find the starting point (first boundary pixel)
    [start_row, start_col] = find(bwperim(binary_image), 1, 'first');
    if isempty(start_row)
        chain_code = [];
        return;
    end

    % Define direction offsets for 8-neighborhood
    directions = [ -1,  0;  % North
                   -1,  1;  % Northeast
                    0,  1;  % East
                    1,  1;  % Southeast
                    1,  0;  % South
                    1, -1;  % Southwest
                    0, -1;  % West
                   -1, -1]; % Northwest

    % Initialize variables
    chain_code = [];
    current_row = start_row;
    current_col = start_col;
    prev_dir = 0; % Start direction

    % Begin tracing
    while true
        found_next = false;
        for i = 0:7
            dir_idx = mod(prev_dir + i, 8) + 1;
            next_row = current_row + directions(dir_idx, 1);
            next_col = current_col + directions(dir_idx, 2);

            % Check bounds
            if next_row < 1 || next_row > size(binary_image, 1) || ...
               next_col < 1 || next_col > size(binary_image, 2)
                continue;
            end

            % Check if the neighbor is part of the perimeter
            if binary_image(next_row, next_col)
                % Add direction to chain code
                chain_code = [chain_code, dir_idx - 1];
                % Move to the next pixel
                current_row = next_row;
                current_col = next_col;
                % Update previous direction
                prev_dir = mod(dir_idx + 4 - 1, 8); % Opposite direction
                found_next = true;
                break;
            end
        end

        % If we returned to the starting point, exit
        if current_row == start_row && current_col == start_col && ~isempty(chain_code)
            break;
        end

        % If no neighbor found, tracing is complete
        if ~found_next
            break;
        end
    end
end

function perimeter = calculate_perimeter_from_chain_code(chain_code)
    % Define distances for horizontal/vertical (1) and diagonal (sqrt(2))
    horizontal_vertical_distance = 1;
    diagonal_distance = sqrt(2);

    % Initialize perimeter
    perimeter = 0;

    % Loop through the chain code and calculate the total distance
    for i = 1:length(chain_code)
        if mod(chain_code(i), 2) == 0
            % Horizontal or vertical move
            perimeter = perimeter + horizontal_vertical_distance;
        else
            % Diagonal move
            perimeter = perimeter + diagonal_distance;
        end
    end
end
