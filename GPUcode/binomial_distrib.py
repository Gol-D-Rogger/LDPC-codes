import math

def binom_pmf(k: int, n: int, p: float)-> float:

Compute P(X k) for X Binomial(n, p) using log-domain arithmetic.

 

if k < 0 or k > n:

return 0.0

if p == 0.0:

return 1.0 if k = 0 else 0.0

if p == 1.0:

return 1.0 if k = n else 0.0

log_pmf =

math. lgamma(n 1)

math. lgamma(k 1)

math. lgamma(n k +1)

k math. log(p)

+ (n- k) * math.log(1.0 - p)

)


return math. exp(log_pmf)

if __name__== "__main__":

total 0

for k_i in range(530, 38912+1): res binom_pmf(p=0.01, k=k_i, n=38272) total + res

print(total)