from .mcelice_bitmap import BitMatrix, BitVector, McElice
from .isd_attack import InformationSetDecoding
import random
from itertools import combinations


class LeeBrickellAttack(InformationSetDecoding):
    """
    Lee-Brickell ISD Variante: Erlaubt p Fehler im Information Set.
    """

    def __init__(self):
        super().__init__()

    def generate_error_vectors(self, length, weight):
        """
        Generiere alle Fehlervektoren der Länge length mit Hamming-Gewicht weight.
        Liefert BitVector-Objekte.
        """
        if weight > length:
            return
        for positions in combinations(range(length), weight):
            error_vec = BitVector(length)
            for pos in positions:
                error_vec.set_bit(pos, 1)
            yield error_vec

    

    def attack(self, ciphertext, t, p, max_attempts=1000000, verbose=True):
        """
        Lee-Brickell ISD Angriff

        Args:
            ciphertext: BitVector - empfangenes verschlüsseltes Wort (Länge n)
            error_weight: int - erwartete Gesamtfehleranzahl t
            p: int - erlaubte Fehler im Information Set (0 = klassisches Prange-ISD)
            max_attempts: int - maximale Anzahl zufälliger Informationssets (Schleifen)
            verbose: bool
        Returns:
            BitVector oder None
        """
        if not isinstance(ciphertext, BitVector):
            raise TypeError("ciphertext muss BitVector sein")

        if self.G_pub is None:
            raise ValueError("G_pub muss gesetzt sein (set_G_pub aufrufen)")

        attempt = 0
        total_checks = 0

        # defensive: p darf nicht größer als k sein
        if p < 0:
            raise ValueError("p muss >= 0 sein")
        if p > self.k:
            p = self.k

        while attempt < max_attempts:
            attempt += 1

            # 1) Zufälliges Informationsset I (k Positionen)
            positions = random.sample(range(self.n), self.k)

            # 2) k×k-Teilmatrix G_I extrahieren
            G_sub = self.extract_submatrix(self.G_pub, positions)  # erwartet: k x k BitMatrix

            # 3) Prüfen, ob invertierbar. Falls nicht: nächstes I
            try:
                G_inv = G_sub.inverse()   # erwartet: k x k BitMatrix (oder wirf Fehler)
            except (ValueError, ZeroDivisionError):
                # nicht invertierbar -> überspringen
                if verbose and (attempt % 10000 == 0):
                    print(f"Versuch {attempt}: G_sub nicht invertierbar (Checks {total_checks})")
                continue

            # 4) y_I extrahieren (k-bit BitVector)
            y_sub = self.extract_subvector(ciphertext, positions)

            # 5) Schleife über i = 0..p
            for i in range(p + 1):
                # Vorfilter: wenn Kombinationsanzahl sehr groß, evtl. skip/limit
                # Generiere alle e_I mit Gewicht i
                for e_I in self.generate_error_vectors(self.k, i):
                    total_checks += 1

                    # (y_I - e_I) in GF(2) ist XOR; angenommen BitVector.add macht XOR
                    y_I_corrected = y_sub.add(e_I)

                    # Baue 1 x k Matrix aus y_I_corrected (für Multiplikation mit G_inv)
                    y_corr_mat = BitMatrix(1, self.k)
                    for idx in range(self.k):
                        y_corr_mat.set_bit(0, idx, y_I_corrected.get_bit(idx))

                    # x_candidate_matrix = (1 x k) * (k x k) = (1 x k)
                    x_candidate_matrix = y_corr_mat.multiply(G_inv)

                    # Erfolgs-Kandidat als BitVector der Länge k
                    x_candidate = BitVector(self.k)
                    for idx in range(self.k):
                        x_candidate.set_bit(idx, x_candidate_matrix.get_bit(0, idx))

                    # 6) globaler Test: wt(y - x'G') == t ?
                    is_correct, wt = self.check_candidate(ciphertext, x_candidate, self.G_pub, t)

                    if is_correct:
                        if verbose:
                            print(f"✓ Erfolg nach {attempt} Informationssets und {total_checks} Checks!")
                            print(f"  Fehlervektorgewicht im Information Set: {i}")
                        return x_candidate

                # optional: wenn viele Kombinationen für dieses i, kann man hier abbrechen/loggen
                if verbose and i % 5 == 0 and i != 0:
                    print(f" Informationsset-Versuch {attempt}: bis p={i} geprüft (Checks {total_checks})")

            if verbose and attempt % 1000 == 0:
                print(f"Versuch {attempt}: {total_checks} Kandidaten geprüft")

        if verbose:
            print(f"✗ Angriff fehlgeschlagen nach {max_attempts} Informationssets ({total_checks} Checks)")
        return None
