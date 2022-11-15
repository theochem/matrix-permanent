import numpy as np

from scipy.special import comb


M_ROWS = 65
N_COLS = 65


a = np.zeros((M_ROWS, N_COLS), dtype=int)

for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        a[i, j] = comb(i, j, exact=True)

a_string = ",\n    ".join(map(str, a.reshape(a.size)))


print(
f"""#ifndef PERMANENT_BINOM_H
#define PERMANENT_BINOM_H


#include <stdint.h>


#define BINOM(N, K) binom_table[{N_COLS} * N + K]


static const int64_t binom_table[] = {{
    {a_string}
}};


#endif /* PERMANENT_BINOM_H */""")
