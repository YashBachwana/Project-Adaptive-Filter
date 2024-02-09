clc ; clear all ; close all ;
input_length = 5000 ;
SNR = 30 ;

input = rand(1,input_length) - 0.5;
% input = 2 * input;

% Calculate the mean of the input
% mean_input = mean(input);

% Adjust the mean to 0
% input = input - mean_input;
system_noise = awgn(input,SNR)-input ;

numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
denominator_coefficients = [1, -1.99, 1.572, -0.4583];
system_output = sin(filter(numerator_coefficients, denominator_coefficients, input)) ;
% system_output = g(input) ;
% for i = 1:input_length
%     system_output(i) = f(input(i)) ;
% end
% system_output = system_output

x = -2:0.2:2 ;
x_out = sin(filter(numerator_coefficients, denominator_coefficients, x)) ;
y = (-2:0.2:2) ;
z = -2:0.2:2 ;
figure(1) ; hold on ; 

plot(x,y,'r') ; hold on ; plot(x,x_out) ; hold on ; 
plot(x, y, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; 
pause(0) ;

Q = 21 ;
P = 2 ; del_x = 0.2 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
m = 1 ; 

for i = 1 : input_length
    cla ;
    i    
    u = (input(i)/del_x) - floor(input(i)/del_x);
    parameter_array_u = [(u)^2,u,1];
    index = floor(input(i)/del_x) + (Q - 1)/2;
    
    out = parameter_array_u * C * y(index : index + 2)' ;

    error = system_output(i) - out ;

    y(index:index + 2) = y(index : index + 2)' + 0.025 * error * C' * parameter_array_u' ;
    % m = m + 0.025 * error * parameter_array_u * C * x(index:index + 2)' ; 
    plot(x,z) ; hold on ; 
    plot(x,y,'r') ; hold on ; plot(x,x_out) ; hold on ; 
    plot(x, y, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ;
    plot(x(index:index + 2), y(index:index + 2), 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 6) ;
    
    pause(0) ;
end 

function physical_output = f(x)
    numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
    denominator_coefficients = [1, -1.99, 1.572, -0.4583];
    physical_output = sin(filter(numerator_coefficients, denominator_coefficients, x)) ;
end 