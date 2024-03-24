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

% system_output = g(input) ;
numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
denominator_coefficients = [1, -1.99, 1.572, -0.4583];
system_output = sin(filter(numerator_coefficients, denominator_coefficients, input)) ; %+ system_noise ;
% for i = 1:input_length
%     system_output(i) = f(input(i)) ;
% end
% system_output = system_output
x_ = -2.5 : 0.05 : 2.5 ;
x_out = sin(filter(numerator_coefficients, denominator_coefficients, x_)) ;

x1 = -2.5:0.05:-1.5 ; x2 = -2:0.05:-1 ; x3 = -1.5:0.05:-0.5 ; x4 = -1:0.05:0 ; x5 = -0.5:0.05:0.5 ;
x6 = 0:0.05:1 ; x7 = 0.5:0.05:1.5 ; x8 = 1:0.05:2 ; x9 = 1.5:0.05:2.5 ;
x = {x1,x2,x3,x4,x5,x6,x7,x8,x9} ;

% y1 = x1 ; y2 = x2 ; y3 = x3 ; y4 = x4 ; y5 = x5 ; y6 = x6 ; y7 = x7 ; y8 = x8 ; y9 = x9 ; y10 = x10 ;
Y = x ; Z = x ;
% z1 = x1 ; z2 = x2 ; z3 = x3 ; z4 = x4 ; z5 = x5 ; z6 = x6 ; z7 = x7 ; z8 = x8 ; z9 = x9 ; z10 = x10 ;
figure(1) ; hold on ; 

for i = 1 : 9
    plot(x{i} , Y{i}) ; hold on ; plot(x_,x_out,'y') ; hold on ; 
    plot(x{i}, Y{i}, 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 4) ; hold on ; 
end 

pause(0) ;

Q = (101:-20:-72)';
P = 2 ; del_x = 0.05 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
mu = 1 ; n = 9 ; 

for i = 1 : input_length
    cla ;
    xlim([-2.5 2.5])
    ylim([-5 5])
    i    
    u = (input(i)/del_x) - floor(input(i)/del_x);
    parameter_array_u = [(u)^2,u,1];
    
    for j = 1 : n                           
        index = floor(input(i)/del_x) + (Q(j) - 1)/2 + 1 ;
        if (index >= 1 && index <= 18)
            out = parameter_array_u * C * Y{j}(index : index + 2)' ;
            error = system_output(i) - out ;
            Y{j}(index:index + 2) = Y{j}(index : index + 2)' + mu * error * C' * parameter_array_u' ;
            plot(x{j},Z{j}) ; hold on ; 
            plot(x{j},Y{j}) ; hold on ; plot(x_,x_out,'y') ; hold on ; 
            plot(x{j}, Y{j}, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;
            plot(x{j}(index:index + 2), Y{j}(index:index + 2), 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 6) ; hold on ; 
        end 
    end 
    pause(0) ;
end 

for j = 1 : n
    plot(x{j},Z{j}) ; hold on ; 
    plot(x{j},Y{j}) ; hold on ; plot(x_,x_out,'y') ; hold on ; 
    plot(x{j}, Y{j}, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;      
end 


% function system_output = g(x)
%     input_length = length(x) ;
%     system_output = zeros(1,input_length) ;
%     for i = 1:input_length
%         system_output(i) = f(x(i)) ;
%     end
% end 
% 
% function physical_output = f(x)
%     numerator_coefficients = [0.0154, 0.0462, 0.0462, 0.0154] ; 
%     denominator_coefficients = [1, -1.99, 1.572, -0.4583];
%     physical_output = sin(filter(numerator_coefficients, denominator_coefficients, x)) ;
% end 