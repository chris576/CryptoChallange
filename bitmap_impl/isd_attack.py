"""
McElice Attack Implementation: Information Set Decoding Attack

Die Idee besteht darin, zufällig k Spalten der Generatormatrix G' zu wählen
und zu hoffen, dass diese Positionen des Geheimtextes keine Fehler enthalten.

Bitmap-Version: Optimiert mit gmpy2 für schnelle Bit-Operationen
"""
from .mcelice_bitmap import (
    McElice, BitMatrix, BitVector, 
    decode_syndrome_bitmap, parity_check_matrix,
    bit_xor, bit_and
)
import random   


class InformationSetDecoding(McElice):

    def __init__(self):
        super().__init__()
        self.G_pub = None
        self.k = 0
        self.n = 0

    def hamming_weight(self, v):
        """Hamming-Gewicht eines BitVectors"""
        if isinstance(v, BitVector):
            return v.weight()
        elif isinstance(v, list):
            return sum(v)
        else:
            raise TypeError("v muss BitVector oder Liste sein")
    
    def calculate_complexity(self, n, k, t):
        """
        Berechne Wahrscheinlichkeit, dass ein zufällig gewähltes Information Set
        keine Fehler enthält.
        
        Args:
            n: Gesamtlänge des Codewortes
            k: Größe des Information Sets
            t: Anzahl der Fehler
        
        Rechnung der durchschnittlichen Laufzeit: 
        = 0(k^3(n über k) / (n-t über k))
        
        Returns:
            float - Wahrscheinlichkeit (0.0 bis 1.0)
        """
        from math import comb
        if t < 0 or t > n:
            raise ValueError("t muss zwischen 0 und n liegen")
        if k > n:
            raise ValueError("k darf nicht größer als n sein")
        if k > n - t:
            print("Warnung: k > n - t, Wahrscheinlichkeit ist 0")
            return 0.0  # Unmöglich, ein fehlerfreies IS zu wählen
        
        numerator = comb(n - t, k)
        denominator = comb(n, k)
        return k**3 * (numerator / denominator)
    
    def set_G_pub(self, G):
        """
        Setze öffentliche Generatormatrix
        
        Args:
            G: BitMatrix - Public Key
        """
        if not isinstance(G, BitMatrix):
            raise TypeError("G muss BitMatrix sein")
        
        self.G_pub = G
        self.k = G.rows
        self.n = G.cols
        return self
    
    def extract_submatrix(self, matrix, positions):
        """
        Extrahiere Teilmatrix aus gegebenen Spaltenpositionen
        
        Args:
            matrix: BitMatrix
            positions: Liste von Spaltenindizes
        
        Returns:
            BitMatrix mit ausgewählten Spalten
        """
        rows = matrix.rows
        cols = len(positions)
        submatrix = BitMatrix(rows, cols)
        
        for i in range(rows):
            for j, col_idx in enumerate(positions):
                submatrix.set_bit(i, j, matrix.get_bit(i, col_idx))
        
        return submatrix
    
    def extract_subvector(self, vector, positions):
        """
        Extrahiere Teilvektor aus gegebenen Positionen
        
        Args:
            vector: BitVector
            positions: Liste von Indizes
        
        Returns:
            BitVector mit ausgewählten Elementen
        """
        length = len(positions)
        subvector = BitVector(length)
        
        for i, pos in enumerate(positions):
            subvector.set_bit(i, vector.get_bit(pos))
        
        return subvector
    
    def matrix_vector_multiply(self, matrix, vector):
        """
        Matrix-Vektor-Multiplikation über GF(2)
        
        Args:
            matrix: BitMatrix (k x n)
            vector: BitVector (k)
        
        Returns:
            BitVector (n)
        """
        if matrix.cols != vector.length:
            raise ValueError(f"Inkompatible Dimensionen: {matrix.cols} != {vector.length}")
        
        result = BitVector(matrix.rows)
        
        for i in range(matrix.rows):
            val = 0
            for j in range(matrix.cols):
                val = bit_xor(val, bit_and(matrix.get_bit(i, j), vector.get_bit(j)))
            result.set_bit(i, val)
        
        return result
    
    def check_candidate(self, cipher, candidate, G_pub, t):
        """
        Prüft, ob ein rekonstruierter Kandidat korrekt ist.

        Args:
            cipher: BitVector - empfangenes Wort (Länge n)
            candidate: BitVector - Kandidat für Klartext (Länge k)
            G_pub: BitMatrix - öffentliche Generatormatrix G' (k x n)
            expected_weight: int - erwartetes Fehlergewicht (meist 0 für fehlerfreie Decodierung)
        
        Returns:
            (bool, int) - ob korrekt, und das Hamming-Gewicht des Fehlers
        """
        # c' = candidate * G_pub (Zeile * Matrix)
        candidate_matrix = BitMatrix(1, candidate.length)
        for j in range(candidate.length):
            candidate_matrix.set_bit(0, j, candidate.get_bit(j))
        
        c_prime_matrix = candidate_matrix.multiply(G_pub)
        c_prime = BitVector(G_pub.cols)
        for j in range(G_pub.cols):
            c_prime.set_bit(j, c_prime_matrix.get_bit(0, j))
        
        # r = cipher XOR c' (Fehlervektor)
        r = cipher.add(c_prime)
        
        # Fehlergewicht prüfen
        wt = r.weight()
        return wt == t, wt
    
    def attack(self, ciphertext, t, max_attempts=1000000, verbose=True):
        """
        Information Set Decoding Angriff
        
        Args:
            ciphertext: BitVector - empfangenes verschlüsseltes Wort
            expected_error_weight: int - erwartetes Fehlergewicht nach Decodierung (Standard: 0)
            num_errors: int - bekannte Anzahl Fehler in der Übertragung (z.B. 11)
                              Wenn > 0, werden bevorzugt fehlerfreie Positionen gewählt
            max_attempts: int - maximale Anzahl Versuche
            verbose: bool - Ausgabe von Fortschritt
        
        Returns:
            BitVector - dekodierte Nachricht oder None
        
        Algorithmus mit Fehlerberücksichtigung:
        1. Wenn num_errors bekannt: Wähle k Positionen so, dass möglichst wenig Fehlerpositionen dabei sind
        2. Extrahiere k×k Teilmatrix aus G_pub an diesen Positionen
        3. Falls invertierbar: Berechne Kandidat x' = G_inv @ y_sub
        4. Prüfe ob Fehlergewicht korrekt ist (sollte = num_errors sein)
        5. Falls ja: Rückgabe, sonst wiederholen
        """
        if not isinstance(ciphertext, BitVector):
            raise TypeError("ciphertext muss BitVector sein")
        
        if self.G_pub is None:
            raise ValueError("G_pub muss gesetzt sein (set_G_pub aufrufen)")
        
        # Wenn Fehleranzahl bekannt: Erhöhe erwartetes Gewicht
        if t > 0:
            if verbose:
                print(f"Erwarte {t} Fehler in der Übertragung")
                self.complexity = self.calculate_complexity(self.n, self.k, t)
                print(f"Geschätzte Anzahl der Versuche: {self.complexity}")                
                
        
        attempt = 0
        import time
        import statistics
        times = []
        while attempt < max_attempts:
            attempt += 1
            # 1. Zufällig k Positionen auswählen
            # Strategie: Wenn wir num_errors kennen, versuchen wir k Positionen zu wählen,
            # die wahrscheinlich keine Fehler enthalten
            # Da wir nicht wissen WO die Fehler sind, bleibt es zufällig, aber wir
            # erwarten, dass bei vielen Versuchen einige Information Sets fehlerfrei sind
            positions = random.sample(range(self.n), self.k)
            
            t1 = time.time()
            
            # 2. Extrahiere k×k Teilmatrix
            G_sub = self.extract_submatrix(self.G_pub, positions)
            
            # 3. Prüfen, ob G_sub invertierbar ist
            try:
                G_inv = G_sub.inverse()
            except (ValueError, ZeroDivisionError):
                # Nicht invertierbar, neue Positionen wählen
                t2 = time.time()
                times.append(t2 - t1)
                if verbose and attempt % 10000 == 0:
                    print(f"Versuch {attempt}: Nicht invertierbar")
                continue
            
            # 4. Extrahiere Teilvektor des Ciphertexts
            y_sub = self.extract_subvector(ciphertext, positions)
            
            # 5. Berechne Kandidat: x' = G_inv @ y_sub
            # y_sub als Matrix (1 x k) behandeln
            y_sub_matrix = BitMatrix(1, self.k)
            for j in range(self.k):
                y_sub_matrix.set_bit(0, j, y_sub.get_bit(j))
            
            x_candidate_matrix = y_sub_matrix.multiply(G_inv)
            x_candidate = BitVector(self.k)
            for j in range(self.k):
                x_candidate.set_bit(j, x_candidate_matrix.get_bit(0, j))
            
            # 6. Fehlervektor berechnen und Hamming-Gewicht prüfen
            is_correct, wt = self.check_candidate(ciphertext, x_candidate, self.G_pub, t)
                        
            t2 = time.time()
            times.append(t2 - t1)
            
            if is_correct:
                if verbose:
                    print(f"✓ Erfolg nach {attempt} Versuchen!")
                    print(f"Median Zeit pro Versuch: {statistics.median(times):.6f} Sekunden")
                    print(f"Geschätzte Zeit laut geschätzten Versuchen * Median-Zeit: {self.complexity * statistics.median(times):.2f} Sekunden")
                    print(f"Tatsächliche Zeit: {sum(times):.2f} Sekunden")
                return x_candidate
        
        if verbose:
            print(f"✗ Angriff fehlgeschlagen nach {max_attempts} Versuchen")
        return None