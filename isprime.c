#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#define CNT(x) (sizeof(x)/sizeof(*x)/!(sizeof(x)%sizeof(*x)))

/* check for square; return sqrt if square else 0 */
uint32_t is_square(uint64_t nn)
{
  assert(nn % 2 != 0);

  static unsigned char inv0[128] = {
    1, 171, 167, 205, 143, 5, 73, 253, 31, 75, 57, 45, 175, 101, 215, 93, 63,
    21, 231, 115, 207, 197, 9, 67, 95, 117, 7, 237, 239, 219, 233, 227, 127,
    213, 217, 77, 241, 123, 55, 125, 159, 203, 71, 83, 209, 27, 169, 35, 191,
    107, 153, 243, 177, 69, 119, 195, 223, 11, 135, 109, 145, 165, 105, 157,
    255, 85, 89, 51, 113, 251, 183, 3, 225, 181, 199, 211, 81, 155, 41, 163,
    193, 235, 25, 141, 49, 59, 247, 189, 161, 139, 249, 19, 17, 37, 23, 29,
    129, 43, 39, 179, 15, 133, 201, 131, 97, 53, 185, 173, 47, 229, 87, 221,
    65, 149, 103, 13, 79, 187, 137, 61, 33, 245, 121, 147, 111, 91, 151, 99
  };

  if (nn % 8 != 1)
    return 0;

  uint32_t inv = inv0[nn / 8 & 127];

  inv *= (3 - inv * inv * (uint32_t)nn) / 2;

  uint64_t r = nn * inv;
  r *= (3 - inv * r) / 2;

  if (r & (1ull << 32))
    r = (uint32_t)(0 - r);
  else
    r = (uint32_t)r;

  if (r * r == nn)
    return (uint32_t)r;

  return 0;

}

/* a * b % m */
uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m)
{
  assert (a < m);
  assert (b < m);
#ifndef __x86_64__
  return (unsigned __int128)a * b % m;
#else
  uint64_t quot, rem;
  asm("mulq %3\n\tdivq %4" : "=a"(quot), "=&d"(rem) : "%a"(a), "r"(b), "r"(m) : "cc");
  return rem;
#endif
}

/* (a + b) % m; no overflow when a + b > 2**64 */
uint64_t addmod(uint64_t a, uint64_t b, uint64_t m)
{
  assert (a < m);
  assert (b < m);
  return (b >= m - a) ? b - (m - a) : b + a;
}

/* 2**power % modul */
uint64_t pow2mod(uint64_t power, uint64_t modul)
{
  uint64_t result;
  uint64_t prb;

  assert(modul > 2);
  assert(power > 0);

  prb = 1ull << (63 - __builtin_clzll(power));
  result = 2;
  while ((prb >>= 1) != 0)
  {
      result = mulmod(result, result, modul);

      if (power & prb)
         result = addmod(result, result, modul);
  }

  return result;
}

/* miller-rabin test for 2 */
int maybeprime(uint64_t n)
{
  assert(n % 2 == 1);
  uint64_t m1 = n - 1;
  unsigned k = __builtin_ctzll(m1);
  uint64_t r = pow2mod(m1 >> k, n);
  if (r == 1 || r == m1)
    return 1;

  while (--k)
  {
     r = mulmod(r, r, n);
     if (r == 1)
       return 0;
     if (r == m1)
       return 1;
  }
  return 0;
}

/* jacobi symbol */
int jacobi(uint64_t a, uint64_t n)
{
    assert(n % 2 == 1);
    assert(a < n);
    int t = 2;
    uint64_t r;

    while (a > 1)
    {
        while (a % 4 == 0)
            a >>= 2;
        if (a % 2 == 0)
        {
            a >>= 1;
            r = n % 8;
            if (r == 3 || r == 5)
                t ^= 2;
        }

        r = n;
        n = a;
        a = r;
        if ((a & 3 & n) == 3)
            t ^= 2;

        a = a % n;
    }

    if (a == 1 || n == 1)
        return t - 1;
    else
        return 0;
}

static unsigned small_prime[] = {
  3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
  61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

/* least quadratic non-residue */
unsigned min_nonres(uint64_t n)
{
  assert(n % 2 == 1);

  /* check for n % 8 in 3, 5 */
//  if ((n ^ (n >> 1)) & 2)
//    return 2;
  assert(n % 8 != 3 && n % 8 != 5);

  uint32_t inv = (n >> 1) & 1;

  if (n % 3 == 2 - inv)
    return 3;
  if (n % 5 == 2 || n % 5 == 3)
    return 5;

  unsigned i;

  // 7 .. 61
  static uint64_t mask64[] = {
     0x17,
     0x23b,
     0x161b,
     0x1a317,
     0x30af3,
     0x5335f,
     0x13d122f3,
     0x121d47b7,
     0x165e211e9bLL,
     0x1b382b50737LL,
     0x35883a3ee53LL,
     0x4351b2753dfLL,
     0x12dd703303aed3LL,
     0x22b62183e7b92bbLL,
     0x1713e6940a59f23bLL
  };

  for (i = 0; i < CNT(mask64); ++i)
  {
     unsigned j = i + 2; //skip 3, 5
     uint8_t mod = (uint8_t)(n % small_prime[j]);
     if (((small_prime[j] >> 1) & inv) != !((1LL << mod) & mask64[i]))
        return small_prime[j];
  }

  // 67 .. 127
  static unsigned __int128 mask128[] = {
     ((unsigned __int128)0x3ULL<<64)|0x59c281ba27ebc653ULL,
     ((unsigned __int128)0x1ULL<<64)|0x164729716b1d977fULL,
     ((unsigned __int128)0x1ebULL<<64)|0x22c742790b8d135fULL,
     ((unsigned __int128)0x130bULL<<64)|0x409e755186fd2f37ULL,
     ((unsigned __int128)0x26873ULL<<64)|0xa80b1372fea31e9bULL,
     ((unsigned __int128)0x1b3c3b9ULL<<64)|0x2a6b59502770f37ULL,
     ((unsigned __int128)0x1eb628b47ULL<<64)|0x6067981b8b451b5fULL,
     ((unsigned __int128)0x1391b7f0d3ULL<<64)|0x552a832c3fb6273ULL,
     ((unsigned __int128)0x16380e9115ULL<<64)|0xbd964257768fe397ULL,
     ((unsigned __int128)0x27816ea9821ULL<<64)|0x633397be6a897e1bULL,
     ((unsigned __int128)0x1752639f4e85ULL<<64)|0xb003685cbe7192bbULL,
     ((unsigned __int128)0x1a75e89ae2121ULL<<64)|0xf33e1211d645eb97ULL,
     ((unsigned __int128)0x172a099c419697f1ULL<<64)|0x7016967dc66fab17ULL
  };

  for (i = 0; i < CNT(mask128); ++i)
  {
     unsigned j = i + 2 + CNT(mask64);
     uint8_t mod = (uint8_t)(n % small_prime[j]);
     if (((small_prime[j] >> 1) & inv) != !(((unsigned __int128)1 << mod) & mask128[i]))
        return small_prime[j];
  }

  if (is_square(n) != 0)
    return 0; //is square

  uint64_t k;
  for (k = 131; ; k += 6)
  {
     if(jacobi(k, n) == -1)
       return k;
     if(jacobi(k + 2, n) == -1)
       return k + 2;
  }
}

uint64_t binary_gcd(uint64_t a, uint64_t b)
{
  /* first must be odd */
  assert (a % 2 == 1);

  b >>= __builtin_ctzll(b);

  /* after right shift high bit will be zero */
  b >>= 1;
  a >>= 1;

  uint64_t c;
  while ((c = a - b) != 0)
  {
     c = (uint64_t)llabs((int64_t)c);
     b = (a + b - c) >> 1; /* min(a, b) */
     a = (c >> 1) >> __builtin_ctzll(c);
  }
  return a + a + 1; /* restore low bit */
}

/* modular powering of (x + y√d) */
uint64_t pp1(uint64_t p, uint64_t m, uint64_t n)
{
  uint64_t d = p * (p - 2);
  uint64_t a = p - 1;
  uint64_t prb = 1ull << (63 - __builtin_clzll(m));
  uint64_t rx = a, ry = 1;

  while ((prb >>= 1) != 0)
  {
      // y' = 2*x*y
      ry = mulmod(rx, ry, n);
      ry = addmod(ry, ry, n);

      // x' = 2*x*x-1
      rx = mulmod(rx, rx, n);
      rx = addmod(rx, rx, n);
      rx = (rx ? rx : n) - 1;

      if (m & prb)
      {
          uint64_t t = addmod(mulmod(a, rx, n), mulmod(ry, d, n), n);
          ry = addmod(mulmod(a, ry, n), rx, n);
          rx = t;
      }
  }
  return rx;
}

uint64_t mul2mod(uint64_t a, uint64_t m)
{
  assert (a < m);
  return (a >= m - a) ? a - (m - a) : a + a;
}

/* modular powering of (x + y√2) */
uint64_t pp2(uint64_t m, uint64_t n)
{
  uint64_t prb = 1ull << (63 - __builtin_clzll(m));
  uint64_t rx = 3, ry = 2;

  while ((prb >>= 1) != 0)
  {
      ry = mul2mod(mulmod(rx, ry, n), n);
      rx = mul2mod(mulmod(rx, rx, n), n);
      rx = (rx ? rx : n) - 1;
      if (m & prb)
      {
          uint64_t t = addmod(addmod(mul2mod(rx, n), rx, n), mul2mod(mul2mod(ry, n), n), n);
          ry = addmod(addmod(mul2mod(ry, n), ry, n), mul2mod(rx, n), n);
          rx = t;
      }
  }
  return rx;
}


/* return a % b */
uint64_t reminder(unsigned __int128 a, uint64_t b)
{
#ifndef __x86_64__
  return a % b;
#else
  uint64_t result;
  asm("divq %3" : "=d"(result): "a"((uint64_t)a), "d"((uint64_t)(a >> 64) % b), "r"(b) : "cc");
  return result;
#endif
}

int isprime(uint64_t n)
{
  if (n == 2)
     return 1;
  if (n % 2 != 1)
     return 0;

  const static __int128 primeprod = (__int128)3*5*7*11*13*17*19*23*29*31*37*41*43*47*53*59*61*67*71*73*79*83*89*97*101;
  uint64_t gcd = binary_gcd(n, reminder(primeprod, n));

  if (gcd != 1)
  {
     if (gcd != n || n > 101)
        return 0;
     unsigned i;
     for (i = 0; small_prime[i] <= n; ++i)
        if (small_prime[i] == n)
           return 1;
     return 0;
  }

  if (!maybeprime(n))
    return 0;

  if (n < 2047) // min base 2 strong pseudoprime
      return 1;


  uint64_t m = (n & 3) == 3 ? (n + 1) / 4 : (n + 1) / 2;
  uint64_t rx;

  /* check for n % 8 in 3, 5 */
  if ((n ^ (n >> 1)) & 2)
     rx = pp2(m, n);
  else
  {
     uint64_t p = min_nonres(n);
     if (p == 0) //is square
        return 0;
     rx = pp1(p, m, n);
  }
  return rx == ((n & 3) == 3 ? 0 : n - 1);
}

