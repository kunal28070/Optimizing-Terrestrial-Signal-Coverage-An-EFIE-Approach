% Define constants
Eamp = 1.0;
Epsilon_0 = 8.854e-12;
Epsilon_d = 2.47912e-11;
Mu_0 = 4 * pi * 1e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0);
GrossStep= 10.0;
f = 150e6;
Lambda = c / f;
DeltaX = Lambda / 4.0;
Omega = 2.0 * pi * f;
j = 1i;
GrossNoSteps = 70;
Beta_0 = Omega * sqrt(Mu_0 * Epsilon_0);
Eta_0 = sqrt(Mu_0 / Epsilon_0);
TOL = 10e-15;
NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);
Xsource = 0.0;
Ysource = 442.0;
I = 1.0;
PI = 3.14159265358979323846;
EXP = 2.718281828;

% Load terrain profile data
data = load('X.04');
X = data(:, 1);
Y = data(:, 2);

% Initialize arrays
ModJ = zeros(1, NoLinesubs);
ModEt = zeros(1, NoLinesubs);

Et = zeros(1, NoLinesubs);
% Sigma = zeros(1, NoLinesubs);

intendedRange = 1500;


x_discrete = 0:DeltaX:intendedRange;
NoLinesubs = length(x_discrete);
y_discrete = interp1(X, Y, x_discrete, 'linear', 'extrap');


E = zeros(NoLinesubs, 1);
Zmat = zeros(NoLinesubs, NoLinesubs);
J = zeros(NoLinesubs, 1);

% Compute E matrix
for j = 1:NoLinesubs
    distance_j = sqrt((x_discrete(j) - 0)^2 + (y_discrete(j) - 442)^2);
    E(j) = besselh(0, 2, Beta_0 * distance_j);
end

% Compute Z matrix 
for i = 1:NoLinesubs
    for j = 1:NoLinesubs
        distance_ij = sqrt((x_discrete(i) - x_discrete(j))^2 + (y_discrete(i) - y_discrete(j))^2);
        if i == j
            Zmat(i, j) = DeltaX * (Beta_0/4) * (1 - (2/i) * log(1.781 * Beta_0 * DeltaX / (4 * e_const)));
        else
            Zmat(i, j) = DeltaX * (Beta_0/4) * besselh(0, 2, Beta_0 * distance_ij);
        end
    end
end

for j = 1:NoLinesubs
    sum_ZJ = 0;
    for i = 1:(j-1)
        sum_ZJ = sum_ZJ + J(i) * Zmat(j, i);
    end
    J(j) = (E(j) - sum_ZJ) / Zmat(j, j);
end

% Plotting
figure;
plot(x_discrete, abs(J));
xlabel('Position (m)');
ylabel('Surface current J');
title('plot of J using the inverse matrix solution at 970 MHz');
grid on;
 
% p

% Calculate path loss and corresponding electric field for each point up to 700 meters

coutput = fopen('E2.dat', 'w');
for index = 1:NoLinesubs
    Et(index) = 0;
    for n = 1:index
        Et(index) = Et(index) + (J(n) * R_p_q(n, n+1, DeltaX, GrossStep, Y) * Z_1(R_surf_obs(n, index, DeltaX, GrossStep, Y), Beta_0, Omega, Epsilon_0));
    end
    Et(index) = EiRad(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource), index, Beta_0, Omega, Epsilon_0) - Et(index);
    fprintf(coutput, '%f  %f\n', DeltaX * (index - 1), 20.0 * log10(abs(Et(index)) / sqrt(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource))));
end
fclose(coutput);
 data2 = load('E2.dat');
  figure; % Create a new figure window
  plot(data2(:,1), data2(:,2));
  xlabel('x coordinates'); % Customize the label
  ylabel('|E(x)|'); % Customize the label
  title('Full wave solution for E(x) at 150 MHz ');
% function definitions start
function result = x(a, DeltaX)    
    result = double(a)* DeltaX;
end

function result = R_source_p(p, DeltaX, GrossStep, Ysource, Xsource, Y)
    result = sqrt(((Xsource-x(p, DeltaX))*(Xsource-x(p, DeltaX)))+((Ysource-y(p, DeltaX,GrossStep, Y))*(Ysource-y(p, DeltaX,GrossStep, Y))));
end

function result = y(a, DeltaX, GrossStep, Y)
    
    Temp=(a*DeltaX)/GrossStep; %tells u what number plate you're on plus a bit
    Index= int64(Temp); %tells u what plate u are on
    Prop=Temp- double(Index); % proportion in x (which is hopefully the same as the proportion in Y!)
    
    s= Y(Index + 1)+(Prop*(Y(Index+2)-Y(Index+1)));
    result = s;
end


function result =  R_source_obs(p, DeltaX, GrossStep, Y, Xsource, Ysource)
    result = (sqrt(((Xsource-x(p, DeltaX))*(Xsource-x(p, DeltaX)))+((Ysource-y(p, DeltaX, GrossStep, Y)-2.4)*(Ysource-y(p, DeltaX, GrossStep, Y)-2.4))));
end

function result = R_p_q(p, q, DeltaX, GrossStep, Y)
    result =(sqrt(((x(q, DeltaX)-x(p, DeltaX))*(x(q, DeltaX)-x(p, DeltaX)))+((y(q, DeltaX, GrossStep, Y)-y(p, DeltaX, GrossStep, Y))*(y(q, DeltaX, GrossStep, Y)-y(p, DeltaX,GrossStep, Y)))));
end

function result = R_surf_obs(p, q, DeltaX, GrossStep, Y)
    result = (sqrt(((x(q, DeltaX)-x(p, DeltaX))*(x(q, DeltaX)-x(p, DeltaX)))+(((y(q, DeltaX, GrossStep, Y)+2.4)-y(p, DeltaX, GrossStep, Y))*((y(q, DeltaX, GrossStep, Y)+2.4)-y(p, DeltaX, GrossStep, Y)))));
end



%complex
function result = EiRad( dist, ~, Beta_0, Omega, Epsilon_0)

    complex E;

    E= -((Beta_0*Beta_0)/(4.0*Omega*Epsilon_0))*H02(Beta_0*dist);


    result = E;
end


function result = H02(Arg)
    complex H;
    H=complex(besselj(0,Arg),-bessely(0,Arg));
    result = (H);
end




function result = Z(p, q, Beta_0, Epsilon_0, DeltaX, GrossStep, Y, Omega)
    complex H;
    
    H=((Beta_0*Beta_0)/(4.0*Omega*Epsilon_0))*H02(Beta_0*R_p_q(p,q,  DeltaX, GrossStep, Y));
    
    
    result = H;
end

function result = Z_1(R, Beta_0, Omega, Epsilon_0)
    complex H;
    H=((Beta_0*Beta_0)/(4.0*Omega*Epsilon_0))*H02(Beta_0*R);
    result = H;
end



% ??
function result = Zself(i, DeltaX, GrossStep,Y, Beta_0, Epsilon_0, Omega, PI, EXP)
double Linesubln;
complex H;

Linesubln=R_p_q(i,i+1, DeltaX, GrossStep, Y);

H=complex(((Beta_0*Beta_0)/(4.0*Omega*Epsilon_0))*Linesubln,-((Beta_0*Beta_0)/(4.0*Omega*Epsilon_0))*((2.0*Linesubln)/PI)*log((1.781*Beta_0*Linesubln)/( 4.0*EXP)));

result = (H);
end




function result = Exp(d)

result = exp(d);
end


function result = cplx_Exp(d)

complex result;

result=complex(cos(d),-sin(d));

end

 