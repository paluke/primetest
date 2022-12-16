import math

def jacobi(a, n):
    assert(n % 2 == 1)
    a %= n
    t = 1
    while a > 1:
        while a & 3 == 0:
            a = a >> 2
        if a & 1 == 0:
            a = a >> 1
            r = n & 7
            if r == 3 or r == 5:
                t = -t
        if a & 3 == 3 and n & 3 == 3:
            t = -t
        a, n = n % a, a
    if a == 1 or n == 1:
        return t
    else:
        return 0

def pairpow(a, m, d, n):
    r = (1, 0)
    while True:
        if m % 2 == 1:
            r = ((r[0] * a[0] + d * r[1] * a[1]) % n, (r[0] * a[1] + a[0] * r[1]) % n)
        m = m >> 1
        if m == 0:
            break
        a = ((2 * a[0] * a[0] + n - 1) % n, (2 * a[0] * a[1]) % n)
    return r

def isprime(n):
    if n == 2:
        return True
    if (n & 1) == 0:
        return False

    if 3 <= (n & 7) <= 5:
        p = 4
        # check 2 is quadratic nonresidue
        if pow(2, n >> 1, n) != n - 1:
            return False
    else:
        if pow(2, n >> 1, n) != 1:
            return False
        m = math.isqrt(n)
        if m * m == n:
            return False
        # seek for least quadratic nonresidue
        p = 1
        j = 1
        while j == 1:
            p = p + 2
            j = jacobi(p, n)
            if j == 0:
                return False
    # final test
    a = pairpow((p - 1, 1), (n + 1) >> 1, p * (p - 2), n)
    if a[0] != n - 1:
        return False
    return True

