import numpy as np

from scipy.special import comb


M_ROWS = 65
N_COLS = 65


binom_table = np.zeros((M_ROWS, N_COLS), dtype=int)

for i in range(binom_table.shape[0]):
    for j in range(binom_table.shape[1]):
        binom_table[i, j] = comb(i, j, exact=True)

binom_table_str = ",\n    ".join(map(str, binom_table.reshape(-1)))


print(
f"""#ifndef PERMANENT_BINOM_H
#define PERMANENT_BINOM_H


#define BINOM(N, K) binom_table[{N_COLS} * N + K]


static const long binom_table[] = {{
    {binom_table_str}
}};


#endif /* PERMANENT_BINOM_H */""")
