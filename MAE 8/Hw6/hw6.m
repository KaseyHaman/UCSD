clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 6;

%%Problem 1
p1a = 0;
for n = 1:40
    p1a = p1a + ((2.^n)/factorial(n));
end

p1b = 0;
for m = 0:40
    for n = 0:40
    p1b = p1b + (1/(3.^(m+n)));
    end
end

p1c = 0;
for m = 0:40
    for n = 0:m
    p1c = p1c + (1/(3.^(m+n)));
    end
end

p1d = 0;
for l = 1:40
    for m = 1:40
        for n = 1:40
          p1d = p1d + 1/((2.^l).*(2.^m).*(2.^n));
        end
    end
end

p1e = 0;
for l = 1:5
    for m = 1:5
        for n = l:m
          p1e = p1e + 1;
        end
    end
end

p1f = 1;
for n = 1:1000
    p1f = p1f.*((4.*(n.^2))/((4.*(n.^2))-1));
end

%%Problem 2
p2a = 0;
while exp(1)/factorial(p2a+1) >= 1e-7
    p2a = p2a + 1;
end

p2b = 0; 
while (2^(p2b+1))*(p2b+1) <= 1e7
    p2b = p2b + 1;
end

%%Problem 3 Double Check

h = 10;
Traveld = zeros(1, 10);
for b = 2:10
    h = h.*(3/4);
    Traveld(b) = Traveld(b-1) + h.*2;
end
Traveld(:) = Traveld(:) + 10;
p3a = 7;

h = 10;
Traveld = 10;
while Traveld < 59.99
    h = h.*(3/4);
    Traveld = Traveld + h.*2;
end
p3b = h;

%%Problem 4 function strfind
load('stringA.mat');
how_counter = 0;
are_counter = 0;
for_counter = 0;
for n = 1:length(stringA)-2
	switch stringA(n:n+2)
		case 'how'
			how_counter = how_counter + 1;
        case 'are'
			are_counter = are_counter + 1;
		case 'for'
			for_counter = for_counter + 1;
  end
end
p4a = how_counter;
p4b = are_counter;
p4c = for_counter;

many_counter = 0;
time_counter = 0;
loop_counter = 0;
for n = 1:length(stringA)-3
	switch stringA(n:n+3)
		case 'many'
			many_counter = many_counter + 1;
        case 'time'
			time_counter = time_counter + 1;
		case 'loop'
			loop_counter = loop_counter + 1;
  end
end
p4d = many_counter;
p4e = time_counter;
p4f = loop_counter;


%%Problem 5

%abs(x-pi)/abs(pi) * 100
cvpos = pi + (10^(-5)); 
cvneg = pi - (10^(-5));

p5a = zeros(1, 200000);
p5a(1) = 1;
k = 200000;
for n = 1:k
    p5a(n+1) = p5a(n) + (((-1).^n)/(2*n+1));
end
p5a(:) = p5a(:) * 4;
a = find(cvneg< p5a & p5a< cvpos);
p5a = p5a(a(1));
p5b = a(1) - 1;

p5c = zeros(1, 100);
p5c(1) = 1;
k = 100;
for n = 1:k
    p5c(n+1) = p5c(n) + (((-3).^(-n))/(2*n+1));
end
p5c(:) = p5c(:) * sqrt(12);
a = find(cvneg< p5c & p5c< cvpos);
p5c = p5c(a(1));
p5d = a(1) - 1 ;

p5e = 'The series in part (c, d) converges faster';



