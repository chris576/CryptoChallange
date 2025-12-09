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
    
    def calculate_success_probability(self, n, k, num_errors):
        """
        Berechne Wahrscheinlichkeit, dass ein zufällig gewähltes Information Set
        keine Fehler enthält.
        
        Args:
            n: Gesamtlänge des Codewortes
            k: Größe des Information Sets
            num_errors: Anzahl der Fehler
        
        Returns:
            float - Wahrscheinlichkeit (0.0 bis 1.0)
        """
        # Wahrscheinlichkeit = C(n-num_errors, k) / C(n, k)
        # = (n-num_errors)! / k! * (n-num_errors-k)! * k! * (n-k)! / n!
        # = Produkt_{i=0}^{k-1} (n-num_errors-i) / (n-i)
        
        if num_errors >= n - k + 1:
            return 0.0  # Unmöglich ein fehlerfreies Information Set zu finden
        
        prob = 1.0
        for i in range(k):
            prob *= (n - num_errors - i) / (n - i)
        
        return prob
    
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
    
    def check_candidate(self, cipher, candidate, G_pub, expected_weight):
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
        return wt == expected_weight, wt
    
    def attack(self, ciphertext, expected_error_weight=0, num_errors=0, max_attempts=1000000, verbose=True):
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
        if num_errors > 0 and expected_error_weight == 0:
            expected_error_weight = num_errors
            if verbose:
                print(f"Erwarte {num_errors} Fehler in der Übertragung")
                success_prob = self.calculate_success_probability(self.n, self.k, num_errors)
                print(f"Wahrscheinlichkeit für fehlerfreies Information Set: {success_prob:.6f}")
                expected_attempts = 1.0 / success_prob if success_prob > 0 else float('inf')
                print(f"Erwartete Anzahl Versuche: {expected_attempts:.0f}\n")
        
        attempt = 0
        while attempt < max_attempts:
            attempt += 1
            
            # 1. Zufällig k Positionen auswählen
            # Strategie: Wenn wir num_errors kennen, versuchen wir k Positionen zu wählen,
            # die wahrscheinlich keine Fehler enthalten
            # Da wir nicht wissen WO die Fehler sind, bleibt es zufällig, aber wir
            # erwarten, dass bei vielen Versuchen einige Information Sets fehlerfrei sind
            positions = random.sample(range(self.n), self.k)
            
            # 2. Extrahiere k×k Teilmatrix
            G_sub = self.extract_submatrix(self.G_pub, positions)
            
            # 3. Prüfen, ob G_sub invertierbar ist
            try:
                G_inv = G_sub.inverse()
            except (ValueError, ZeroDivisionError):
                # Nicht invertierbar, neue Positionen wählen
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
            is_correct, wt = self.check_candidate(ciphertext, x_candidate, self.G_pub, expected_error_weight)
            
            if verbose and (attempt % 10000 == 0 or is_correct):
                print(f"Versuch {attempt}: Fehlergewicht {wt} (erwartet: {expected_error_weight})")
            
            if is_correct:
                if verbose:
                    print(f"✓ Erfolg nach {attempt} Versuchen!")
                return x_candidate
        
        if verbose:
            print(f"✗ Angriff fehlgeschlagen nach {max_attempts} Versuchen")
        return None