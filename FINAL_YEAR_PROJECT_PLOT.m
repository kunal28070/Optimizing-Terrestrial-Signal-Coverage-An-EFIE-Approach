% Define constants
Eamp = 1.0;
Epsilon_0 = 8.854e-12;
Epsilon_d = 2.47912e-11;
Mu_0 = 12.56637061e-7;
c = 1.0 / sqrt(Mu_0 * Epsilon_0);
GrossStep= 10.0;
f = 970e6;
Lambda = c / f;
DeltaX = Lambda / 4.0;
Omega = 2.0 * pi * f;
j = 1i;
GrossNoSteps = 150;
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
J = zeros(1, NoLinesubs);
Et = zeros(1, NoLinesubs);
% Sigma = zeros(1, NoLinesubs);




% Forward scattering
J(1) = EiRad(R_source_p(1, DeltaX, GrossStep, Ysource, Xsource, Y), 1, Beta_0, Omega, Epsilon_0) / Zself(1, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP);

for p = 1:NoLinesubs
    SUM = 0;
    for q = 1:p-1
        SUM = SUM + R_p_q(q, q+1, DeltaX, GrossStep, Y) * Z(p, q, Beta_0, Epsilon_0, DeltaX, GrossStep, Y, Omega) * J(q);
    end
    J(p) = (EiRad(R_source_p(p, DeltaX, GrossStep, Ysource, Xsource, Y), p, Beta_0, Omega, Epsilon_0) - SUM) / Zself(p, DeltaX, GrossStep, Y, Beta_0, Epsilon_0, Omega, PI, EXP);
end



% Write current values to file
%dlmwrite('J.dat', [DeltaX*(0:NoLinesubs-1); abs(J)]');
writematrix([DeltaX*(0:NoLinesubs-1); abs(J)]', 'J.dat')
data1 = load('J.dat');
  figure; % Create a new figure window
  plot(data1(:,1), data1(:,2));
  xlabel('x coordinates'); % Customize the label
  ylabel('J(x)'); % Customize the label
  title('Plot of J(x)  at 970 MHz');
  
% Calculate total electric field above the surface
coutput = fopen('E.dat', 'w');
for index = 1:NoLinesubs
    Et(index) = 0;
    for n = 1:index
        Et(index) = Et(index) + (J(n) * R_p_q(n, n+1, DeltaX, GrossStep, Y) * Z_1(R_surf_obs(n, index, DeltaX, GrossStep, Y), Beta_0, Omega, Epsilon_0));
    end
    Et(index) = EiRad(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource), index, Beta_0, Omega, Epsilon_0) - Et(index);
    fprintf(coutput, '%f  %f\n', DeltaX * (index - 1), 20.0 * log10(abs(Et(index)) / sqrt(R_source_obs(index, DeltaX, GrossStep, Y, Xsource, Ysource))));
end
fclose(coutput);
 data2 = load('E.dat');
  figure; % Create a new figure window
  plot(data2(:,1), data2(:,2));
  xlabel('x coordinates'); % Customize the label
  ylabel('|E(x)|'); % Customize the label
  title('Plot of E(x) at 970 MHz ');

% Constants for path loss calculation
h_receiver = 2; % Receiver antenna height in meters
C_suburban = calculate_correction_factor_suburban(h_receiver, f);
data2 = load('E.dat');

% Ensure x_discrete and data2 have the same length
x_discrete = data2(:, 1).';
Et = data2(:, 2).';








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

 