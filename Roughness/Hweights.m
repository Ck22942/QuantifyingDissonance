H2 = [	0	0
        17    0.8
        23	0.95
        25	0.975
        32	1
        37	0.975
        48	0.9
        67    0.8
        90	0.7
        114   0.6
        171   0.4
        206	0.3
        247   0.2
        294	0.1
        358	0 ];

x = H2(: , 1);
y = H2(:, 2);

function y = my_function(x, a, b, c, d)

    y = a .* exp(-b .* x) - c .* exp(-d.*(x));
end

n = length(x);

x =  0 : 1 : n;
a = 1; b = 3.5; c = 1; d = 5.0; % Set parameters
y = my_function(x, a, b, c, d); % Evaluate the function
plot(x, y);
xlabel('x'); ylabel('f(x)');
title('Function Plot');
grid on;

