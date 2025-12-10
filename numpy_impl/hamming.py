import numpy as np
from typing import Sequence
from galois import GF2 

# Generatormatrix G
# Fehlerhaftes Codeword R
# Ursprüngliches Codeword mittels Bruteforce finden.
# Bruteforce über alle möglichen Informationsvektoren x, so dass xG dem empfangenen Wort R am nächsten kommt.
# Rückgabe des rekonstruierten Codewords C = xG sowie des Informationsvektors x.
def decode_bruteforce(G: Sequence[Sequence[int]], R: Sequence[int]) -> tuple[np.ndarray | None, np.ndarray | None]:
    G = GF2(G)
    R = GF2(R).flatten()

    if G.ndim != 2:
        raise ValueError("Generatormatrix G muss zweidimensional sein.")

    k, n = G.shape

    if R.size != n:
        raise ValueError(f"Codeword R muss Länge {n} haben, hat aber {R.size}.")

    for i in range(2**k):
        info_bits = GF2.Zeros(k)
        bits = list(np.binary_repr(i, width=k))
        for idx, bit in enumerate(bits):
            info_bits[idx] = int(bit)

        codeword = info_bits @ G % 2
        error_vector = (R + codeword) % 2

        if int(error_vector.sum()) <= (n - k) // 2:
            return np.array(codeword, dtype=int), np.array(info_bits, dtype=int)

    return None, None

def ParityCheckMatrix(G: np.ndarray) -> np.ndarray:
    """Berechnet die Parity-Check-Matrix H für eine Generatormatrix G über GF(2)."""
    G_bin = np.array(G, dtype=int) % 2
    return _gf2_null_space(G_bin)


def _gf2_null_space(matrix: np.ndarray) -> np.ndarray:
    """Berechnet eine Basis des Nullraums über GF(2) mittels Gauß-Elimination."""
    A = GF2(matrix)
    m, n = A.shape
    A = A.copy()

    pivot_columns: list[int] = []
    row = 0
    for col in range(n):
        if row >= m:
            break
        pivot_candidates = np.where(A[row:, col] == 1)[0]
        if pivot_candidates.size == 0:
            continue
        pivot = pivot_candidates[0] + row
        if pivot != row:
            A[[row, pivot]] = A[[pivot, row]]
        for r in range(m):
            if r != row and A[r, col] == 1:
                A[r, :] ^= A[row, :]
        pivot_columns.append(col)
        row += 1

    free_columns = [c for c in range(n) if c not in pivot_columns]
    if not free_columns:
        return np.zeros((0, n), dtype=int)

    basis = []
    for free_col in free_columns:
        vector = GF2.Zeros(n)
        vector[free_col] = 1
        for r, pivot_col in enumerate(pivot_columns):
            if A[r, free_col] == 1:
                vector[pivot_col] = 1
        basis.append(np.array(vector, dtype=int))

    return np.array(basis, dtype=int)


def decode_syndrome(H, codeword):
    """Dekodiert ein empfangenes Wort mittels Syndrom-Dekodierung über GF(2)."""
    H = np.array(GF2(H), dtype=int)
    r = np.array(GF2(codeword), dtype=int).flatten()

    if H.ndim != 2:
        raise ValueError("Parity-Check-Matrix H muss zweidimensional sein.")

    n = H.shape[1]
    if r.size != n:
        raise ValueError(f"Empfangenes Wort muss Länge {n} haben, hat aber {r.size}.")

    k = n - H.shape[0]

    # 1. Berechne s = Hr^T
    s = (H @ r) % 2
    # 2. Ist s = 0, so decodiere y := r = nn(r)
    if np.all(s == 0):
        return r, r[:k]
    # 3. Prüfe für alle Fehlervektoren mit e \in F^n mit wt(e) = 1, ob s = He^T. Falls ja, dekodiere 
    # y:= r +e = nn(r)
    for i in range(n):
        e = np.zeros(n, dtype=int)
        e[i] = 1
        syndrome_e = np.dot(H, e) % 2
        if np.array_equal(s, syndrome_e):
            y = (r + e) % 2
            return y, y[:k]
        
    # 4. Fehe wie in 3. vor für alle Vektoren e mit wt(e) = 2
    for i in range(n):
        for j in range(i + 1, n):
            e = np.zeros(n, dtype=int)
            e[i] = 1
            e[j] = 1
            syndrome_e = np.dot(H, e) % 2
            if np.array_equal(s, syndrome_e):
                y = (r + e) % 2
                return y, y[:k]
    # 5. Ist stets s != He^T oder die Dekodierung liefert kein Codewort, 
    # so sind mehr als [d-1/2] Fehler aufgetreten: Return too many transmission errors
    return None, None
    