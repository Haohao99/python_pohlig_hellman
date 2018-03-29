import timeit

def showDiscreteLog(p, a):
    print('p: ' + str(p))
    print('a: ' + str(a))
    print()
    x = 1
    while x < p:
        print('x: ' + str(x) + '\tb: ' + str(pow(a, x, p)))
        x = x + 1

def CRT(dict_mod_num):
    M = 1
    for num, mod in dict_mod_num.items():
        M = M * mod
    ans = 0
    for num, mod in dict_mod_num.items():
        m = (int(M / mod)) % mod
        ans = ans + (num * inverse(m, mod) * (int(M / mod)))
    return ans % M

def inverse(a, n):
    t = 0
    newt = 1
    r = n
    newr = a
    while newr != 0:
        quotient = int(r / newr)
        t, newt = newt, t - quotient * newt
        r, newr = newr, r - quotient * newr
    if t < 0:
        t = t + n
    return t

def findDiscreteLog(a, B, p):
    primes = prime_factors(p - 1)
    CRT_dict = {}
    for prime, power in primes.items():
        exp = 1
        x_list = []
        B_sub = B
        while exp <= power:
            q = int((p - 1) / pow(prime, exp))
            Bpow = pow(B_sub, q, p)
            apow = pow(a, int((p - 1) / prime), p)
            k = 0
            found = False
            while not found:
                if pow(apow, k, p) == Bpow:
                    x_list.append(k)
                    # Note that this only checks up to q - 1
                    found = True
                if not found:
                    k = k + 1
            B_sub = (B_sub * pow(inverse(a, p), k * pow(prime, exp - 1))) % p
            exp = exp + 1
        exp = 0
        ans_mod_prime = 0
        while exp < power:
            ans_mod_prime = ans_mod_prime + (pow(prime, exp) * x_list[exp])
            exp = exp + 1
        CRT_dict[ans_mod_prime % pow(prime, power)] = pow(prime, power)
    return CRT(CRT_dict)

def prime_factors(n):
    i = 2
    factors = {}
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            try:
                factors[i] = factors[i] + 1
            except KeyError:
                factors[i] = 1
    if n > 1:
        try:
            factors[n] = factors[n] + 1
        except KeyError:
            factors[n] = 1
    return factors

if __name__ == '__main__':
    #showDiscreteLog(13, 7)
    #dict_num_mod = {65 : 99, 2 : 98, 51 : 97, 10 : 95}
    #print(CRT(dict_num_mod))
    #print(prime_factors(2 * 2 * 2 * 3 * 5 * 7 *7 * 7* 17))
    #start = timeit.default_timer()
    print('1: '+ str(findDiscreteLog(3, 8576, 53047))) # problem 1
    #stop = timeit.default_timer()
    #print('Time Elapsed: '+ str(stop - start) + ' s')
    print('2: '+ str(findDiscreteLog(3, 24, 31))) # problem 2
    print('3a: '+ str(findDiscreteLog(2, 3925, 3989))) # problem 3
    print('3a: '+ str(findDiscreteLog(2, 1046, 3989))) # problem 3
    print('3b: '+ str(findDiscreteLog(2, (3925 * 1046) % 3989, 3989))) # pr 3
    print('4: '+ str(findDiscreteLog(11, 2, 1201))) # problem 4



