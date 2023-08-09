clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 1;

%%Problem 1
p1a = pi/0;
p1b = 0/0;
p1c = sqrt(-4*pi);
p1d = cosd(75);
p1e = sin (pi/3);
p1f = (1234)^5;
p1g = nthroot(512, 9);
p1h = log2(16384);
p1i = log10(1000000);
p1j = log(exp(1));
p1k = atand(1);
p1l = sinh(6);
p1m = atanh(1);

%%Problem 2
p2a = char(32*pi);                              
p2b = char(16);
p2c = double(p2b);               
p2d = int16('Z');
p2e = int32('Z');
p2f = single('Y');
p2g = double('X');
p2h = class(p2a);
p2i = double('Y')*double('Z');
p2j = (double('Y') == int64('Y'));

%%Problem 3
p3a = ('y'=='Y');
p3b = ('y'>'X');
p3c = ('z'<'x');
p3d = (log2(1024) == 10);
p3e = (sin(100*pi)~=0);
p3f = ((3\9 + 1) < 3);
p3g = ((3/9+1) < 3);
a=3; b=4; c=5;
p3h = (a>b && a>c);
p3i = (a<b && a>c);
p3j = (a>b || a>c);
p3k = (a>b || a<c);
p3l = xor(a<b, a<c);

%%Problem 4
p4a = (fix(2.5) == floor(2.5));
p4b = (fix(2.4) == fix(-2.4));
p4c = (fix(2.2) == floor(2.2));
p4d = (fix(-2.2) == floor(-2.2));
p4e = (fix(-2.2) == ceil(-2.2));
p4f = (rem(5,3) == mod(5,3));
p4g = (rem(5, -3) == mod(5, -3));

%%Problem 5
p5a = '+';
p5b = 'rational';
p5c = 'bank';
p5d = 'shortg';

%%Problem 6
p6a = 'upper';
p6b = abs('A'-'a');
a = char(65:90);
b = char(97:122);
c = a-b;
p6c = abs(c) == 32;

%%Problem 7
pounds = 1000; kilos = pounds * 2.2;
p7a = kilos;
fn = 5.6; dynes = fn * 10^5;
p7b = dynes;
ftemp = 212; ctemp = ((ftemp - 32) * 5/9);
p7c = ctemp;
mph = 65; kmph = mph * 1.6093;
p7d = kmph;

%%Problem 8
p8a = linspace(2, 998, 499);
p8b = 1:2:999;
p8c = [p8b p8a];
p8d = length(p8c);
p8e = find(p8c == 500);
p8f = [0 p8c];
p8g = p8f(((end/4) + 1):(end/2));
p8h = p8f(((3*end/4) + 1):end);
p8i = -1:-2:-1999;
p8j = (p8i).^2;
p8k = sum(p8i);
p8l = prod(p8i(:, (end-4):end));
p8m = cumsum(p8i);

%%Problem 9
A = [1 2; 3 4];
B = zeros(2);
C = [5 5; 5 5];
D = fliplr(A);

p9a = [A B C B D; B A C D B; C C A C C; B D C A B; D B C B A];
p9b = sum(p9a(:, 5));
p9c = sum(diag(p9a)) + sum(diag(fliplr(p9a)));
p9d = sum(p9a, "all");
p9e = sum(p9a > 2, "all");



