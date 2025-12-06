import numpy as np
import galois

GF2 = galois.GF(2)

# Lösen eines LGS mit Gauß-Prozedur auf Basis von galois (GF(2))
# Utility Funktion, die ich mit Hilfe von Claude Sonnet 4.5 erstellt habe, 
# weil ich sagemath nicht installiert bekommen habe in meinem Setup und 
# sympy mir zu viel war. Denke numpy reicht vollkommen aus.
def loese_binaeres_lgs(Ab, return_details=False):
    """
    Löst ein lineares Gleichungssystem über GF(2) mit Hilfe von galois.

    Args:
        Ab: Augmentierte Matrix [A|b] als Array oder verschachtelte Listen mit Einträgen 0/1.
        return_details: Wenn True, werden U, b' und Metadaten zurückgegeben.

    Returns:
        - Wenn konsistent und return_details=False:
            (spezielle_loesung, nullraum)
        - Wenn inkonsistent und return_details=False:
            (None, [])
        - Wenn return_details=True:
            {
              'status': 'consistent' | 'inconsistent',
              'U': U,
              'b_trans': bT,
              'pivot_cols': pivot_cols,
              'solution': spezielle_loesung oder None,
              'nullspace': nullraum (Liste) oder []
            }
    """
    Ab_gf = GF2(Ab)
    _, cols = Ab_gf.shape
    n_vars = cols - 1

    Ab_rref = gauss_vorwaerts_elimination(Ab_gf)
    pivot_cols = ermittle_pivotspalten(Ab_rref[:, :n_vars])
    inconsistent = pruefe_inkonsistenz(Ab_rref, n_vars)

    U = np.array(Ab_rref[:, :n_vars], dtype=int)
    b_trans = np.array(Ab_rref[:, -1], dtype=int)

    if inconsistent:
        if return_details:
            return {
                'status': 'inconsistent',
                'U': U,
                'b_trans': b_trans,
                'pivot_cols': pivot_cols,
                'solution': None,
                'nullspace': []
            }
        return None, []

    spezielle_loesung = berechne_spezielle_loesung(Ab_rref, pivot_cols, n_vars)
    nullraum = berechne_nullraum(Ab_rref[:, :n_vars])

    if return_details:
        return {
            'status': 'consistent',
            'U': U,
            'b_trans': b_trans,
            'pivot_cols': pivot_cols,
            'solution': spezielle_loesung,
            'nullspace': nullraum
        }

    return spezielle_loesung, nullraum


def gauss_vorwaerts_elimination(Ab):
    """Führt die Gaußsche Vorwärtselimination über GF(2) auf einer Matrix aus."""
    Ab_gf = GF2(Ab)
    return Ab_gf.row_reduce()


def ermittle_pivotspalten(matrix):
    """Bestimmt Pivotspalten anhand der RREF-Matrix."""
    matrix_int = np.array(GF2(matrix), dtype=int)
    pivot_cols = []
    for row in matrix_int:
        ones = np.nonzero(row)[0]
        if ones.size > 0:
            pivot_cols.append(int(ones[0]))
    return pivot_cols


def pruefe_inkonsistenz(Ab_rref, n_vars):
    """Prüft auf Zeilen der Form [0 ... 0 | 1] und damit Inkonsistenzen."""
    for row in np.array(GF2(Ab_rref), dtype=int):
        if np.all(row[:n_vars] == 0) and row[-1] == 1:
            return True
    return False


def berechne_spezielle_loesung(Ab_rref, pivot_cols, n_vars):
    """Extrahiert eine spezielle Lösung aus der RREF-Matrix."""
    solution = GF2.Zeros(n_vars)
    for row_idx, pivot_col in enumerate(pivot_cols):
        solution[pivot_col] = Ab_rref[row_idx, -1]
    return np.array(solution, dtype=int)


def berechne_nullraum(A):
    """Erzeugt Basisvektoren des Nullraums über GF(2)."""
    nullspace = GF2(A).null_space()
    if nullspace.size == 0:
        return []
    return [np.array(vec, dtype=int) for vec in nullspace]


def invert(matrix):
    gf_matrix = GF2(matrix)
    if gf_matrix.shape[0] != gf_matrix.shape[1]:
        raise ValueError("Matrix muss quadratisch sein, um invertiert zu werden.")

    n = gf_matrix.shape[0]
    augmented = np.hstack((gf_matrix, GF2.Identity(n)))
    rref = augmented.row_reduce()
    linke_haelfte = rref[:, :n]

    if not np.array_equal(np.array(linke_haelfte, dtype=int), np.eye(n, dtype=int)):
        rang = gf_matrix.row_space().shape[0]
        raise np.linalg.LinAlgError(
            f"Matrix ist singulär: Rang {rang} < {n}."
        )

    inverse = rref[:, n:]
    return np.array(inverse, dtype=int)
    
def matrix_vector_mult(A, v):
    gf_A = GF2(A)
    gf_v = GF2(v)
    return np.array(gf_A @ gf_v, dtype=int)


def format_gf2_arrays(arrays):
    """Formatiert GF(2)-Arrays als eine Zeile mit komma-getrennten Einträgen."""
    def format_vector(vector):
        flat = np.array(GF2(vector), dtype=int).ravel().tolist()
        return "[" + ",".join(str(bit) for bit in flat) + "]"

    data = arrays
    if isinstance(data, (list, tuple)) and not isinstance(data, np.ndarray):
        formatted = [format_vector(vec) for vec in data]
        return ",".join(formatted)

    arr = np.array(GF2(data), dtype=int)
    if arr.ndim == 1:
        return format_vector(arr)
    return ",".join(format_vector(row) for row in arr)