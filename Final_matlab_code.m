clear all
% Parameters

% Use these values as a reference
% k1 = 5; k2 = 5;   % Spring constants (N/m)
% m1 = 2; m2 = 2;   % Masses (kg)
% L1 = 2; L2 = 2;   % Unstretched spring lengths (m)
% w1 = 5; w2 = 5;   % Widths of the masses (m)
% h = 0.01;         %step size 


k1 = input('Please enter spring constant 1:\n');
k2 = input('Please enter spring constant 2:\n');
m1 = input('Please enter mass 1:\n');
m2 = input('Please enter mass 2:\n');
L1 = input('Please enter the unstretched length of spring 1:\n');
L2 = input('Please enter the unstretched length of spring 2:\n');
w1 = input('Please enter the width of mass 1:\n');
w2 = input('Please enter the width of mass 2:\n');
h = input('Please enter the step size:\n');


% Initial conditions
i=1;
x1(i) = L1;  
x2(i) = L1 + w1 + L2 + 6; 
v1(i)= 0;  
v2(i) = 0;  
t0 = 0; 
tf = 20; 
t0 = t0:h:tf;  
n_steps = length(t0);

% Initialize arrays for results
x1_euler(i) =0;
x2_euler(i)=0;
v1_euler(i) = 0;
v2_euler(i) =0;

x1_heun(i) = 0;
x2_heun(i) =0;
v1_heun(i) = 0;
v2_heun(i) = 0;

x1_rk4(i) = 0;
x2_rk4(i)= 0;
v1_rk4(i) =0;
v2_rk4(i) =0;

% Set initial values
x1_euler(1) = x1(i); 
x2_euler(1) = x2(i); 
v1_euler(1) = v1(i);
v2_euler(1) = v2(i);

x1_heun(1) = x1(i);
x2_heun(1) = x2(i);
v1_heun(1) = v1(i);
v2_heun(1) = v2(i);

x1_rk4(1) = x1(i);
x2_rk4(1) = x2(i);
v1_rk4(1) = v1(i);
v2_rk4(1) = v2(i);

% Function definitions for derivatives
dx1dt = @(v1) v1;
dx2dt = @(v2) v2;
dv1dt = @(x1, x2) (-k1/m1)*(x1 - L1) + (k2/m1)*(x2 - x1 - w1 - L2);
dv2dt = @(x1, x2) (-k2/m2)*(x2 - x1 - w1 - L2);

% Euler's Method
for i = 1:n_steps-1
    x1_euler(i+1) = x1_euler(i) + h * dx1dt(v1_euler(i));
    x2_euler(i+1) = x2_euler(i) + h * dx2dt(v2_euler(i));
    v1_euler(i+1) = v1_euler(i) + h * dv1dt(x1_euler(i), x2_euler(i));
    v2_euler(i+1) = v2_euler(i) + h * dv2dt(x1_euler(i), x2_euler(i));
    Ea_x1_euler(i)=abs((x1_euler(i+1)-x1_euler(i))/x1_euler(i+1)*100);
    Ea_x2_euler(i)=abs((x2_euler(i+1)-x2_euler(i))/x2_euler(i+1)*100);
     Ea_v1_euler(i)=abs((v1_euler(i+1)-v1_euler(i))/v1_euler(i+1)*100);
     Ea_v2_euler(i)=abs((v2_euler(i+1)-v2_euler(i))/v2_euler(i+1)*100);


end

% Heun's Method
for i = 1:n_steps-1
    % Predictor step
    x1_pred = x1_heun(i) + h * dx1dt(v1_heun(i));
    x2_pred = x2_heun(i) + h * dx2dt(v2_heun(i));
    v1_pred = v1_heun(i) + h * dv1dt(x1_heun(i), x2_heun(i));
    v2_pred = v2_heun(i) + h * dv2dt(x1_heun(i), x2_heun(i));
    
    % Corrector step
    x1_heun(i+1) = x1_heun(i) + (h/2) * (dx1dt(v1_heun(i)) + dx1dt(v1_pred));
    x2_heun(i+1) = x2_heun(i) + (h/2) * (dx2dt(v2_heun(i)) + dx2dt(v2_pred));
    v1_heun(i+1) = v1_heun(i) + (h/2) * (dv1dt(x1_heun(i), x2_heun(i)) + dv1dt(x1_pred, x2_pred));
    v2_heun(i+1) = v2_heun(i) + (h/2) * (dv2dt(x1_heun(i), x2_heun(i)) + dv2dt(x1_pred, x2_pred));
    Ea_x1_heun(i)=abs((x1_heun(i+1)-x1_heun(i))/x1_heun(i+1)*100);
    Ea_x2_heun(i)=abs((x2_heun(i+1)-x2_heun(i))/x2_heun(i+1)*100);
    Ea_v1_heun(i)=abs((v1_heun(i+1)-v1_heun(i))/v1_heun(i+1)*100);
    Ea_v2_heun(i)=abs((v2_heun(i+1)-v2_heun(i))/v2_heun(i+1)*100);
end

% Runge-Kutta Method (4th order)
for i = 1:n_steps-1
   
    k1x1 =  h * dx1dt(v1_rk4(i));
    k1x2 =  h * dx2dt(v2_rk4(i));
    k1v1 =  h * dv1dt(x1_rk4(i), x2_rk4(i));
    k1v2 =  h * dv2dt(x1_rk4(i), x2_rk4(i));
    
    k2x1 = h * dx1dt(v1_rk4(i) + k1v1/2);
    k2x2 = h * dx2dt(v2_rk4(i) + k1v2/2);
    k2v1 = h * dv1dt(x1_rk4(i) + k1x1/2, x2_rk4(i) + k1x2/2);
    k2v2 = h * dv2dt(x1_rk4(i) + k1x1/2, x2_rk4(i) + k1x2/2);
    
    k3x1 = h * dx1dt(v1_rk4(i) + k2v1/2);
    k3x2 = h * dx2dt(v2_rk4(i) + k2v2/2);
    k3v1 = h * dv1dt(x1_rk4(i) + k2x1/2, x2_rk4(i) + k2x2/2);
    k3v2 = h * dv2dt(x1_rk4(i) + k2x1/2, x2_rk4(i) + k2x2/2);
    
    k4x1 = h * dx1dt(v1_rk4(i) + k3v1);
    k4x2 = h * dx2dt(v2_rk4(i) + k3v2);
    k4v1 = h * dv1dt(x1_rk4(i) + k3x1, x2_rk4(i) + k3x2);
    k4v2 = h * dv2dt(x1_rk4(i) + k3x1, x2_rk4(i) + k3x2);
    
    % Update values
    x1_rk4(i+1) = x1_rk4(i) + (k1x1 + 2*k2x1 + 2*k3x1 + k4x1)/6;
    x2_rk4(i+1) = x2_rk4(i) + (k1x2 + 2*k2x2 + 2*k3x2 + k4x2)/6;
    v1_rk4(i+1) = v1_rk4(i) + (k1v1 + 2*k2v1 + 2*k3v1 + k4v1)/6;
    v2_rk4(i+1) = v2_rk4(i) + (k1v2 + 2*k2v2 + 2*k3v2 + k4v2)/6;
    Ea_x1_rk4(i)=abs((x1_rk4(i+1)-x1_rk4(i))/x1_rk4(i+1)*100);
    Ea_x2_rk4(i)=abs((x2_rk4(i+1)-x2_rk4(i))/x2_rk4(i+1)*100);
    Ea_v1_rk4(i)=abs((v1_rk4(i+1)-v1_rk4(i))/v1_rk4(i+1)*100);
    Ea_v2_rk4(i)=abs((v2_rk4(i+1)-v2_rk4(i))/v2_rk4(i+1)*100);
    



end

% Display results
disp('Comparison of Methods: Final Displacements and Velocities');
disp('-------------------------------------------------------');
disp('Method      | x1 (m)   | x2 (m)   | v1 (m/s) | v2 (m/s)');
disp('-------------------------------------------------------');
fprintf('Euler       | %.4f   | %.4f   | %.4f   | %.4f\n', ...
    x1_euler(end), x2_euler(end), v1_euler(end), v2_euler(end));
fprintf('Heun        | %.4f   | %.4f   | %.4f   | %.4f\n', ...
    x1_heun(end), x2_heun(end), v1_heun(end), v2_heun(end));
fprintf('Runge-Kutta | %.4f   | %.4f   | %.4f   | %.4f\n', ...
    x1_rk4(end), x2_rk4(end), v1_rk4(end), v2_rk4(end));
disp('-------------------------------------------------------');
% Display the approximate errors directly for each method
disp('----------------------------------------------------------');
disp('Approximate Errors for Each Method (at Final Time Step)');
disp('----------------------------------------------------------');
disp('Method      | Ea(x1)  | Ea(x2)  | Ea(v1)  | Ea(v2)');
disp('----------------------------------------------------------');

% Euler Method errors
fprintf('Euler       | %.4f%% | %.4f%% | %.4f%% | %.4f%%\n', ...
    Ea_x1_euler(end), Ea_x2_euler(end), Ea_v1_euler(end), Ea_v2_euler(end));

% Heun Method errors
fprintf('Heun        | %.4f%% | %.4f%% | %.4f%% | %.4f%%\n', ...
    Ea_x1_heun(end), Ea_x2_heun(end), Ea_v1_heun(end), Ea_v2_heun(end));

% Runge-Kutta Method errors
fprintf('Runge-Kutta | %.4f%% | %.4f%% | %.4f%% | %.4f%%\n', ...
    Ea_x1_rk4(end), Ea_x2_rk4(end), Ea_v1_rk4(end), Ea_v2_rk4(end));

disp('----------------------------------------------------------');


% Plot results
% Heun's graph overlaps with Runge-Kutta due to very close results
figure;
subplot(2, 1, 1);
plot(t0, x1_euler, 'r', t0, x1_heun, 'g', t0, x1_rk4, 'b');
xlabel('Time (s)'); ylabel('x1 (m)');
legend('Euler', 'Heun', 'RK4'); 
title('Displacement of Mass 1');
grid on;

subplot(2, 1, 2);
plot(t0, x2_euler, 'r', t0, x2_heun, 'g', t0, x2_rk4, 'b');
xlabel('Time (s)'); ylabel('x2 (m)');
legend('Euler', 'Heun', 'RK4'); 
title('Displacement of Mass 2');
grid on;

figure;
subplot(2, 1, 1);
plot(t0, v1_euler, 'r', t0, v1_heun, 'g', t0, v1_rk4, 'b');
xlabel('Time (s)'); ylabel('v1 (m/s)');
legend('Euler', 'Heun', 'RK4'); 
title('Velocity of Mass 1');
grid on;

subplot(2, 1, 2);
plot(t0, v2_euler, 'r', t0, v2_heun, 'g', t0, v2_rk4, 'b');
xlabel('Time (s)'); ylabel('v2 (m/s)');
legend('Euler', 'Heun', 'RK4'); 
title('Velocity of Mass 2');
grid on;
