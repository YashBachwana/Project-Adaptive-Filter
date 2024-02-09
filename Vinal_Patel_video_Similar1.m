clc ; clear all ; close all ;
input_length = 5000 ;
SNR = 30 ;

input = rand(1,input_length) - 0.5;
input = 2 * input;

% Calculate the mean of the input
mean_input = mean(input);

% Adjust the mean to 0
input = input - mean_input;
system_noise = awgn(input,SNR)-input ;

system_output = g(input) ;
% for i = 1:input_length
%     system_output(i) = f(input(i)) ;
% end
% system_output = system_output

x = -2:0.2:2 ;
y = (-2:0.2:2) ;
z = -2:0.2:2 ;
figure(1) ; hold on ; 

plot(x,y,'r') ; hold on ; plot(x,g(x)) ; hold on ; 
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
    
    out = parameter_array_u * C * y(index : index + 2)' + system_noise(i);

    error = f(input(i)) - out ;

    y(index:index + 2) = y(index : index + 2)' + 1 * error * C' * parameter_array_u' ;
    % m = m + 0.025 * error * parameter_array_u * C * x(index:index + 2)' ; 
    plot(x,z) ; hold on ; 
    plot(x,y,'r') ; hold on ; plot(x,g(x)) ; hold on ; 
    plot(x, y, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ;
    plot(x(index:index + 2), y(index:index + 2), 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 6) ;
    
    pause(0) ;
end 

function system_output = g(x)
    input_length = length(x) ;
    system_output = zeros(1,input_length) ;
    for i = 1:input_length
        system_output(i) = f(x(i)) ;
    end
end 

function physical_output = f(x)
    k1 = 40 ; k2 = 30 ; h1 = 5 ; h2 = 0.5 ; c1 = -0.8 ; c2 = 0.5 ;
    physical_output = k1 * exp(-1 * ((x - c1) ^ 2) / (2 * h2^2)) +  k2 * exp(-1 * ((x - c2) ^ 2) / (2 * h1^2));
end 