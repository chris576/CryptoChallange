# McElice Implementierung mit Bitmaps (optimiert mit gmpy2)
# Grundidee: Verwende Python Integers als Bitvektoren für GF(2) Operationen

import gmpy2


def popcount(x):
    """Zähle die Anzahl der gesetzten Bits (Hamming-Gewicht)"""
    return gmpy2.popcount(x)


def bit_and(a, b):
    """Bitweises AND"""
    return a & b


def bit_or(a, b):
    """Bitweises OR"""
    return a | b


def bit_xor(a, b):
    """Bitweises XOR"""
    return a ^ b


def bit_not(a, width):
    """Bitweises NOT (invertiere alle Bits bis width)"""
    mask = (1 << width) - 1
    return a ^ mask


def test_bit(x, pos):
    """Teste ob Bit an Position pos gesetzt ist"""
    return gmpy2.bit_test(x, pos)


def set_bit(x, pos):
    """Setze Bit an Position pos"""
    return int(gmpy2.bit_set(x, pos))


def clear_bit(x, pos):
    """Lösche Bit an Position pos"""
    return int(gmpy2.bit_clear(x, pos))


class BitMatrix:
    """Matrix über GF(2) als Liste von Integers (jedes Bit = ein Matrixelement)"""
    
    def __init__(self, rows, cols, data=None):
        """
        Args:
            rows: Anzahl Zeilen
            cols: Anzahl Spalten
            data: Liste von Integers (jeder Integer = eine Zeile als Bitmap)
                  oder Liste von Listen mit 0/1 Werten
        """
        self.rows = rows
        self.cols = cols
        
        if data is None:
            self.data = [0] * rows
        elif isinstance(data[0], int):
            self.data = data[:]
        else:
            # Konvertiere Liste von Listen zu Integers
            self.data = []
            for row in data:
                row_int = 0
                for i, bit in enumerate(row):
                    if bit:
                        row_int |= (1 << (self.cols - 1 - i))
                self.data.append(row_int)
    
    def get_bit(self, row, col):
        """Hole Bit an Position (row, col)"""
        pos = self.cols - 1 - col
        return int(gmpy2.bit_test(self.data[row], pos))
    
    def set_bit(self, row, col, value):
        """Setze Bit an Position (row, col)"""
        pos = self.cols - 1 - col
        if value:
            self.data[row] = int(gmpy2.bit_set(self.data[row], pos))
        else:
            self.data[row] = int(gmpy2.bit_clear(self.data[row], pos))
    
    def multiply(self, other):
        """Matrix-Multiplikation über GF(2): self @ other"""
        if isinstance(other, BitVector):
            # Matrix * Vektor
            result = BitVector(self.rows)
            for i in range(self.rows):
                # XOR aller Bits wo sowohl Matrix-Zeile als auch Vektor 1 haben
                overlap = self.data[i] & other.data
                result.data[i] = gmpy2.popcount(overlap) & 1
            return result
        elif isinstance(other, BitMatrix):
            # Matrix * Matrix
            result = BitMatrix(self.rows, other.cols)
            for i in range(self.rows):
                for j in range(other.cols):
                    val = 0
                    for k in range(self.cols):
                        val ^= self.get_bit(i, k) & other.get_bit(k, j)
                    result.set_bit(i, j, val)
            return result
        else:
            raise TypeError("Kann nur mit BitVector oder BitMatrix multiplizieren")
    
    def transpose(self):
        """Transponiere Matrix"""
        result = BitMatrix(self.cols, self.rows)
        for i in range(self.rows):
            for j in range(self.cols):
                result.set_bit(j, i, self.get_bit(i, j))
        return result
    
    def inverse(self):
        """Berechne Inverse über GF(2) mit Gauß-Jordan"""
        if self.rows != self.cols:
            raise ValueError("Matrix muss quadratisch sein")
        
        n = self.rows
        # Erstelle erweiterte Matrix [A | I]
        augmented = BitMatrix(n, 2 * n)
        for i in range(n):
            for j in range(n):
                augmented.set_bit(i, j, self.get_bit(i, j))
            augmented.set_bit(i, n + i, 1)  # Identitätsmatrix
        
        # Gauß-Jordan Elimination
        for col in range(n):
            # Finde Pivot
            pivot_row = None
            for row in range(col, n):
                if augmented.get_bit(row, col) == 1:
                    pivot_row = row
                    break
            
            if pivot_row is None:
                raise ValueError("Matrix ist nicht invertierbar")
            
            # Tausche Zeilen
            if pivot_row != col:
                augmented.data[col], augmented.data[pivot_row] = \
                    augmented.data[pivot_row], augmented.data[col]
            
            # Eliminiere alle anderen Einträge in dieser Spalte
            for row in range(n):
                if row != col and augmented.get_bit(row, col) == 1:
                    augmented.data[row] ^= augmented.data[col]
        
        # Extrahiere inverse Matrix (rechte Hälfte)
        inverse = BitMatrix(n, n)
        for i in range(n):
            for j in range(n):
                inverse.set_bit(i, j, augmented.get_bit(i, n + j))
        
        return inverse
    
    def to_list(self):
        """Konvertiere zu Liste von Listen"""
        result = []
        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                row.append(self.get_bit(i, j))
            result.append(row)
        return result
    
    def __repr__(self):
        return f"BitMatrix({self.rows}x{self.cols})"


class BitVector:
    """Vektor über GF(2) als Integer (jedes Bit = ein Vektorelement)"""
    
    def __init__(self, length, data=None):
        """
        Args:
            length: Länge des Vektors
            data: Integer (als Bitmap) oder Liste von 0/1 Werten
        """
        self.length = length
        
        if data is None:
            self.data = [0] * length
        elif isinstance(data, int):
            # Einzelner Integer - konvertiere zu Liste
            self.data = []
            for i in range(length):
                self.data.append((data >> (length - 1 - i)) & 1)
        else:
            self.data = [int(b) for b in data][:length]
    
    def get_bit(self, pos):
        """Hole Bit an Position pos"""
        return self.data[pos]
    
    def set_bit(self, pos, value):
        """Setze Bit an Position pos"""
        self.data[pos] = int(bool(value))
    
    def add(self, other):
        """Addition über GF(2) = XOR"""
        result = BitVector(self.length)
        for i in range(self.length):
            result.data[i] = self.data[i] ^ other.data[i]
        return result
    
    def weight(self):
        """Hamming-Gewicht (Anzahl der 1en)"""
        return sum(self.data)
    
    def to_list(self):
        """Konvertiere zu Liste"""
        return self.data[:]
    
    def __repr__(self):
        return f"BitVector({self.length}, {self.data})"


def parity_check_matrix(G):
    """Berechne Parity-Check-Matrix H für Generatormatrix G über GF(2)"""
    # Vereinfachte Version: Für systematische Codes G = [I_k | P]
    # ist H = [-P^T | I_{n-k}] = [P^T | I_{n-k}] (da -1 = 1 in GF(2))
    
    k = G.rows
    n = G.cols
    m = n - k
    
    # Extrahiere P (die letzten n-k Spalten von G)
    P = BitMatrix(k, m)
    for i in range(k):
        for j in range(m):
            P.set_bit(i, j, G.get_bit(i, k + j))
    
    # H = [P^T | I_m]
    H = BitMatrix(m, n)
    P_T = P.transpose()
    
    for i in range(m):
        for j in range(k):
            H.set_bit(i, j, P_T.get_bit(i, j))
        H.set_bit(i, k + i, 1)  # Identität
    
    return H


def decode_syndrome_bitmap(H, codeword):
    """
    Dekodiert ein empfangenes Wort mittels Syndrom-Dekodierung über GF(2)
    
    Args:
        H: BitMatrix - Parity-Check-Matrix
        codeword: BitVector - empfangenes Wort
    
    Returns:
        (decoded_codeword, message) oder (None, None) bei zu vielen Fehlern
    """
    n = H.cols
    m = H.rows
    k = n - m
    
    # 1. Berechne Syndrom s = H @ r^T
    r_T = BitMatrix(n, 1)
    for i in range(n):
        r_T.set_bit(i, 0, codeword.get_bit(i))
    
    s_matrix = H.multiply(r_T)
    s = BitVector(m)
    for i in range(m):
        s.set_bit(i, s_matrix.get_bit(i, 0))
    
    # 2. Ist s = 0, keine Fehler
    if s.weight() == 0:
        message = BitVector(k)
        for i in range(k):
            message.set_bit(i, codeword.get_bit(i))
        return codeword, message
    
    # 3. Prüfe Fehlervektoren mit Gewicht 1
    for i in range(n):
        e = BitVector(n)
        e.set_bit(i, 1)
        
        # Berechne H @ e^T
        e_T = BitMatrix(n, 1)
        for j in range(n):
            e_T.set_bit(j, 0, e.get_bit(j))
        
        syndrome_e_matrix = H.multiply(e_T)
        syndrome_e = BitVector(m)
        for j in range(m):
            syndrome_e.set_bit(j, syndrome_e_matrix.get_bit(j, 0))
        
        # Vergleiche Syndrome
        if all(s.get_bit(j) == syndrome_e.get_bit(j) for j in range(m)):
            # Korrigiere Fehler
            corrected = codeword.add(e)
            message = BitVector(k)
            for j in range(k):
                message.set_bit(j, corrected.get_bit(j))
            return corrected, message
    
    # 4. Prüfe Fehlervektoren mit Gewicht 2
    for i in range(n):
        for j in range(i + 1, n):
            e = BitVector(n)
            e.set_bit(i, 1)
            e.set_bit(j, 1)
            
            e_T = BitMatrix(n, 1)
            for idx in range(n):
                e_T.set_bit(idx, 0, e.get_bit(idx))
            
            syndrome_e_matrix = H.multiply(e_T)
            syndrome_e = BitVector(m)
            for idx in range(m):
                syndrome_e.set_bit(idx, syndrome_e_matrix.get_bit(idx, 0))
            
            if all(s.get_bit(idx) == syndrome_e.get_bit(idx) for idx in range(m)):
                corrected = codeword.add(e)
                message = BitVector(k)
                for idx in range(k):
                    message.set_bit(idx, corrected.get_bit(idx))
                return corrected, message
    
    # Zu viele Fehler
    return None, None


class McElice:
    """
    McElice Kryptosystem mit Bitmap-Implementierung (ohne numpy/galois)
    
    Attributes:
        G_priv: BitMatrix - Private Generatormatrix
        S: BitMatrix - Private invertierbare Matrix
        P: BitMatrix - Private Permutationsmatrix
        G_pub: BitMatrix - Public Key Generatormatrix
        H: BitMatrix - Parity-Check-Matrix
    """
    
    def __init__(self):
        self.G_priv = None
        self.S = None
        self.P = None
        self.G_pub = None
        self.H = None
    
    def set_G_priv(self, G):
        """
        Setze private Generatormatrix
        
        Args:
            G: Liste von Listen (Matrix über GF(2))
        """
        self.G_priv = BitMatrix(len(G), len(G[0]), G)
        return self
    
    def set_S(self, S):
        """
        Setze private invertierbare Matrix S
        
        Args:
            S: Liste von Listen (quadratische Matrix über GF(2))
        """
        self.S = BitMatrix(len(S), len(S[0]), S)
        return self
    
    def set_P(self, P):
        """
        Setze Permutationsmatrix P
        
        Args:
            P: Liste von Listen (Permutationsmatrix über GF(2))
        """
        self.P = BitMatrix(len(P), len(P[0]), P)
        return self
    
    def generate_public_key(self):
        """
        Berechne Public Key: G_pub = S @ G_priv @ P
        """
        if self.S is None or self.G_priv is None or self.P is None:
            raise ValueError("S, G_priv und P müssen gesetzt sein")
        
        temp = self.S.multiply(self.G_priv)
        self.G_pub = temp.multiply(self.P)
        return self
    
    def generate_H(self):
        """Berechne Parity-Check-Matrix H aus G_priv"""
        if self.G_priv is None:
            raise ValueError("G_priv muss gesetzt sein")
        
        self.H = parity_check_matrix(self.G_priv)
        return self
    
    def encrypt(self, message, error_vector):
        """
        Verschlüssele Nachricht
        
        Args:
            message: Liste von Bits (Klartext)
            error_vector: Liste von Bits (Fehlervektor)
        
        Returns:
            BitVector - verschlüsseltes Codewort
        """
        if self.G_pub is None:
            raise ValueError("G_pub muss berechnet sein (generate_public_key)")
        
        msg_vec = BitVector(len(message), message)
        err_vec = BitVector(len(error_vector), error_vector)
        
        # Berechne message @ G_pub
        msg_matrix = BitMatrix(1, len(message), [message])
        codeword_matrix = msg_matrix.multiply(self.G_pub)
        
        # Konvertiere zu Vektor
        codeword = BitVector(self.G_pub.cols)
        for i in range(self.G_pub.cols):
            codeword.set_bit(i, codeword_matrix.get_bit(0, i))
        
        # Addiere Fehlervektor
        return codeword.add(err_vec)
    
    def decrypt(self, codeword):
        """
        Entschlüssele Codewort
        
        Args:
            codeword: BitVector oder Liste von Bits
        
        Returns:
            Liste von Bits (Klartext) oder None bei Fehler
        """
        if isinstance(codeword, list):
            codeword = BitVector(len(codeword), codeword)
        
        # Schritt 1: y' = y @ P^-1
        P_inv = self.P.inverse()
        
        # Konvertiere codeword zu Matrix
        y_matrix = BitMatrix(1, codeword.length)
        for i in range(codeword.length):
            y_matrix.set_bit(0, i, codeword.get_bit(i))
        
        y_prime_matrix = y_matrix.multiply(P_inv)
        y_prime = BitVector(codeword.length)
        for i in range(codeword.length):
            y_prime.set_bit(i, y_prime_matrix.get_bit(0, i))
        
        # Schritt 2: Dekodiere y' zu c1
        decoded_result = decode_syndrome_bitmap(self.H, y_prime)
        
        if decoded_result[0] is None:
            raise ValueError("Dekodierung fehlgeschlagen: Zu viele Fehler")
        
        c1, x1 = decoded_result
        
        # Schritt 3: x1 ist bereits die Nachricht (aus decode_syndrome)
        # Schritt 4: Berechne x = x1 @ S^-1
        S_inv = self.S.inverse()
        
        x1_matrix = BitMatrix(1, x1.length)
        for i in range(x1.length):
            x1_matrix.set_bit(0, i, x1.get_bit(i))
        
        plaintext_matrix = x1_matrix.multiply(S_inv)
        
        plaintext = []
        for i in range(plaintext_matrix.cols):
            plaintext.append(plaintext_matrix.get_bit(0, i))
        
        return plaintext