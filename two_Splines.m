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
x = -1.5:0.1:1.5 ;
x1 = -1.5:0.1:0.5 ; x2 = -0.5:0.1:1.5 ;
y1 = -1.5:0.1:0.5 ; y2 = -0.5:0.1:1.5 ;
z1 = -1.5:0.1:0.5 ; z2 = -0.5:0.1:1.5 ;
figure(1) ; hold on ; 

plot(x1,y1,'r') ; hold on ; plot(x,g(x),'y') ; hold on ; 
plot(x1, y1, 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 4) ; hold on ; 

plot(x2,y2,'g') ; hold on ; plot(x,g(x),'y') ; hold on ; 
plot(x2, y2, 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 4) ; hold on ; 
pause(0) ;

Q = [31,11] ;
P = 2 ; del_x = 0.1 ;
C = [0.5,-1,0.5;-1,1,0;0.5,0.5,0] ;
mu = 0.025 ;

for i = 1 : input_length
    cla ;
    xlim([-1.5 1.5])
    i    
    u = (input(i)/del_x) - floor(input(i)/del_x);
    parameter_array_u = [(u)^2,u,1];
    index1 = floor(input(i)/del_x) + (Q(1) - 1)/2 + 1;
    if (index1 >= 1 && index1 <= 18)
        a = index1 ;
        out = parameter_array_u * C * y1(index1 : index1 + 2)' ;
    
        error = f(input(i)) - out ;
    
        y1(index1:index1 + 2) = y1(index1 : index1 + 2)' + mu * error * C' * parameter_array_u' ;
        % m = m + 0.025 * error * parameter_array_u * C * x(index:index + 2)' ; 
        plot(x1,z1) ; hold on ; 
        plot(x1,y1,'r') ; hold on ; plot(x,g(x),'y') ; hold on ; 
        plot(x1, y1, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;
        plot(x1(index1:index1 + 2), y1(index1:index1 + 2), 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 6) ;  hold on ;
    end 
    
    index2 = floor(input(i)/del_x) + (Q(2) - 1)/2 + 1;
    if (index2 >= 1 && index2 <= 18)
        b = index2 ;
        out = parameter_array_u * C * y2(index2 : index2 + 2)' ;
    
        error = f(input(i)) - out ;
    
        y2(index2:index2 + 2) = y2(index2 : index2 + 2)' + mu * error * C' * parameter_array_u' ;
        % m = m + 0.025 * error * parameter_array_u * C * x(index:index + 2)' ; 
        plot(x2,z2) ; hold on ; 
        plot(x2,y2,'g') ; hold on ; plot(x,g(x),'y') ; hold on ; 
        plot(x2, y2, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;
        plot(x2(index2:index2 + 2), y2(index2:index2 + 2), 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 6) ; 
    end 

    pause(0) ;
end 

plot(x1,z1) ; hold on ; 
plot(x1,y1,'r') ; hold on ; plot(x,g(x),'y') ; hold on ; 
plot(x1, y1, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;


plot(x2,z2) ; hold on ; 
plot(x2,y2,'g') ; hold on ; plot(x,g(x),'y') ; hold on ; 
plot(x2, y2, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4) ; hold on ;

figure ;
x = -1.5:0.1:1.5 ;
z = g(x) ;
y = zeros(1,30) ; j = 6 ;
for i = 1 : 31
    if (i <= 15)
        y(i) = y1(i) ;
    else 
        y(i) = y2(j) ;
        j = j + 1 ;
    end 

end 
plot(x,z) ; hold on ; 
plot(x,y) ;

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