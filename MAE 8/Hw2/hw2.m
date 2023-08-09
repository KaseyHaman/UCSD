clear all;
close all;
clc;
format long;
name = 'Kasey haman';
id = 'A16978114';
hw_num = 2;

%%Problem 1
p1a = randi([5, 5], 5, 5);
p1b = p1a;
p1b(:, 2)=0;
p1c = p1b';
p1d = rot90(rot90(rot90(p1b)));
p1e = p1c == p1d;

%%Problem 2
A = [1 0; 0 2]; B = zeros(2); C = fliplr(A); D = [3 3; 4 4];
E = [0 5; 6 0]; F = fliplr(E);
p2a = [A B C A B C; B D B B D B; E B F E B F; A B C A B C; B D B B D B; E B F E B F];
p2b = [A B; B D];
p2c = [B C; D B];
p2d = [D B; B F];
p2e = p2a(4:9, 4:9);
p2f = p2a;
p2f([2, 4, 6, 9:11], :)=[];
p2g = [p2a(1), p2a(1, 12); p2a(length(p2a)), p2a(end)];

%%Problem 3
x = 1:10;
y = 10:10:100;
p3a = x.*y;
p3b = x.^(log10(y));
p3c = (sin(y.^x))./(exp(1).^(y./x));
p3d = (x+((exp(1)).^((-y).^x)))./(y+log(x.^(y)));

%%Problem 4
tmp1 = ones(1,100);
n=1:100;
tmp1(n) = (-1).^(n);
tmp2 = 1:2:199;
p4 = sum(tmp1./tmp2);

%%Problem 5
A = [1:3; 4:6; 7:9]; B = [7 8; 9 10]; C = [1 3 11; 2 7 4];
D = [9 4; 8 5; 3 2];
p5a = A^2;
p5b = 'error';
p5c = 'error';
p5d = A*D;
p5e = 'error';
p5f = D*B;
p5g = B*C;
p5h = 'error';
p5i = 'error';
p5j = isequal(C, D);

%%Problem 6
A = [9 8 7; 6 9 7; 1 7 4];
p6a = A;
b = [1; 3; 5];
p6b = b;
x = A\b;
p6c = x;
x = inv(A)*b; 
p6d = x;
p6e = (p6c == p6d);
p6f = (p6c - p6d);
