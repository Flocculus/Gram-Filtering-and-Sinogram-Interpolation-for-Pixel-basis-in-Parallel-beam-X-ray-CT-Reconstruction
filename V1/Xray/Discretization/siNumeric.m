function y = siNumeric_old(xIn)
% based on NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING pg. 258
% computes the integral from -inf to x of sin(t)/t dt
% equivalent to but faster than built-in sinint

MAXITER = 100;
FPMIN = 1e-30;
EPS = eps;

signs = sign(xIn);
xIn = abs(xIn);

crossover = 2;
y = zeros(size(xIn));

%% for large x
cond = xIn>crossover;
x = xIn(cond);



b = 1 + 1i*x;
c = 1/FPMIN;

d = 1 ./ b;
h = d;


% i = 2:MAXITER;
% a = -(i-1)^2;
% b = b + 2*(i-1);

for i = 2:MAXITER
	a = -(i-1)^2;
	b = b + 2;
	
	d = 1 ./ (b + a * d);
	c = b + a ./ c;
	
	del = c .* d;
	h = h .* del;
	
	if all (abs(real(del)-1) + abs(imag(del)) < EPS)
		break
	end
end

h = h .* (cos(x) - 1i * sin(x));
y(cond) = pi/2 + imag(h);
%% for small x
cond = xIn >= sqrt(FPMIN) & xIn <= crossover;
x = xIn(cond);

sum = zeros(size(x));
sums = sum;
sumc = sum;

signCur = 1;
fact = ones(size(x));
odd = true;

for i = 1:MAXITER
	fact = fact .* x/i;
	term = fact / i;
	sum = sum + signCur * term;
	err = term ./ abs(sum);
	if odd
		signCur = -signCur;
		sums = sum;
		sum = sumc;
	else
		sumc = sum;
		sum = sums;
	end
	
	if (err < EPS) break; end
	odd = ~odd;
end

y(cond) = sums;

%% tiny x
y(xIn < sqrt(FPMIN)) = xIn((xIn < sqrt(FPMIN)) );

%% fix signs
y = y .* signs;