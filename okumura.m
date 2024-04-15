% Load Terrain Data
terrain_data = load("X.04");
X = terrain_data(:, 1)';
Y = terrain_data(:, 2)';

% Load data from E1.dat
E1_data = load('E1.dat');
E1_distances = E1_data(:, 1);
E1_values = E1_data(:, 2);

% Constants
c = 3e8; % Speed of light
f = 970e6 / 1e6; % Frequency in MHz
lambda = c / f;
DeltaX = lambda / 4;
intendedRange = 700; % Range in meters

% Define antenna heights (example values)
h_transmitter = 442; % Transmitter antenna height in meters
h_receiver = 2; % Receiver antenna height in meters

% Okumura-Hata model parameters for suburban areas
A = 69.55 + 26.16*log10(f) - 13.82*log10(h_transmitter);
B = 44.9 - 6.55*log10(h_transmitter);
C = calculate_correction_factor(h_receiver, f);

% Calculate electric field loss for each point up to 700 meters
electric_field_loss = zeros(1, 700);
for i = 1:700
    R = i / 1000; % Distance to the point in kilometers
    electric_field_loss(i) = (A + B*log10(R) - C)-70; % Correction factor for suburban area
end

% Create a figure and plot both data sets
figure;

% Plot E1 data
plot(E1_distances, E1_values, 'b', 'DisplayName', 'Electric field from EFIE');
hold on;

% Plot Electric Field Loss and invert its values for visual y-axis inversion
plot(X(X <= 700), -electric_field_loss(X <= 700), 'r', 'DisplayName', 'Electric Field Loss from Okumura-Hata');

% Add labels and title
xlabel('Distance (m)');
ylabel('Metric Value');
title('Comparison of Electric Field Loss and E1 Data');
legend show; % Show legend to identify the plots
grid on;
hold off;

% Correction factor calculation function for suburban areas
function C = calculate_correction_factor(h_receiver, f)
    % h_receiver: Receiver antenna height in meters
    % f: Frequency in MHz
    if f >= 150 && f <= 1500 % Frequency range for which the model is applicable
        C = -2*(log10(f/28))^2 - 5.4 + (1.1*log10(f) - 0.7)*h_receiver - (1.56*log10(f) - 0.8); % Suburban area correction
    else
        error('Frequency out of range for Okumura-Hata model');
    end
end
