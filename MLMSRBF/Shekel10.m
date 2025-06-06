function y = Shekel10(x)

m = 10;

a = [4, 4, 4, 4;
     1, 1, 1, 1;
     8, 8, 8, 8;
     6, 6, 6, 6;
     3, 7, 3, 7;
     2, 9, 2, 9;
     5, 5, 3, 3;
     8, 1, 8, 1;
     6, 2, 6, 2;
     7, 3.6, 7, 3.6];
 
c = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5];
 
tempsum = 0;
for i = 1:m
    denom = (norm(x-a(i,:)))^2 + c(i);
    tempsum = tempsum + (1/denom);
end

y = -tempsum;