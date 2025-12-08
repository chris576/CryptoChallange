# McElice Implementierung mit Hamming Codes
# Grundidee: Starte mit einem leicht zu decodierendem Code (Secret Key) und transformiere diesen iin einen zufälligen Code (Public Key). 

import numpy as np 
import galois
from hamming import decode_syndrome, GF2, ParityCheckMatrix


class McElice():
    """
    MCElice Kryptosystem basierend auf Hamming Codes.

    Attributes:

    """
    def __init__(self):
        pass
    
    def P(self, P) -> 'McElice': 
        self.P = GF2(P)
        return self
    
    def G_Priv(self, G) -> 'McElice':
        """
        Secret Key Generator Matrix
        Restrictions: 
        - n = 2^m
        - d = 2*t+1
        - k = n-mt
        - S ist eine invertrierbare k*k Matrix
        - P ist eine n*n Permutationsmatrix
        - t ist die Fehlerkorrekturkapazität des Codes
        """
        self.G_Priv = GF2(G)
        return self
    
    def G_Pub(self, S, P) -> 'McElice':
        """
        Public Key Generator Matrix
        G_Pub = S * G_Priv * P
        """
        self.S = GF2(S)
        self.P = GF2(P)
        self.G_Pub = self.S @ self.G_Priv @ self.P
        return self

    def S(self, S) -> 'McElice': 
        self.S = GF2(S)
        return self
    
    def H(self, H) -> 'McElice':
        self.H = GF2(H)
        return self

    def Generate_H(self) -> 'McElice': 
        self.H = ParityCheckMatrix(self.G_Priv)
        return self

    def encrypt(self, message, error_vector):
        """
        Klartext message wird verschlüsselt zu: 
        codeword = message * G + error_vector
        """
        msg_gf2 = GF2(message)
        err_gf2 = GF2(error_vector)
        return msg_gf2 @ self.G_Pub + err_gf2

    def decrypt(self, codeword):
        """
        Entschlüsselung des Codeworts:
        1. Berechne: y' = y * P^-1 = xSG + eP^-1 = c1 + e1, c1 in C
        2. Dekodiere y' zu c1 
        3. Löse x1G = c1 mittels Gausscher Elimination
        4. Berechne den Klartext x = x1*S^-1
        """
        # Konvertiere codeword zu GF(2)
        codeword_gf2 = GF2(codeword)
        
        # Schritt 1: y' = y * P^-1
        P_inv = np.linalg.inv(self.P)
        y_prime = codeword_gf2 @ P_inv

        # Schritt 2: Dekodiere y' zu c1
        decoded_result = decode_syndrome(self.H, y_prime)
        
        if decoded_result[0] is None:
            raise ValueError("Dekodierung fehlgeschlagen: Zu viele Übertragungsfehler")
        
        # c1 ist das korrigierte Codewort (ohne Fehler)
        c1 = GF2(decoded_result[0])

        # Schritt 3: Löse x1*G = c1 mittels Gausscher Elimination über GF(2)
        # Das System ist x1 * G = c1, wobei x1 ein Zeilenvektor ist
        # Wir lösen stattdessen G^T * x1^T = c1^T
        k = self.G_Priv.shape[0]
        n = self.G_Priv.shape[1]
        
        # Erstelle erweiterte Matrix [G^T | c1^T]
        G_T = self.G_Priv.T
        c1_column = c1.reshape(-1, 1)
        augmented_matrix = GF2(np.hstack((G_T, c1_column)))
        
        # Gauß-Elimination über GF(2) zu RREF
        num_rows, num_cols = augmented_matrix.shape
        current_row = 0
        
        for col in range(k):  # Für jede Spalte der Matrix G^T
            if current_row >= num_rows:
                break
                
            # Suche Pivot in dieser Spalte
            pivot_found = False
            for row in range(current_row, num_rows):
                if augmented_matrix[row, col] == 1:
                    # Tausche Zeilen
                    if row != current_row:
                        augmented_matrix[[current_row, row]] = augmented_matrix[[row, current_row]]
                    pivot_found = True
                    break
            
            if not pivot_found:
                continue
            
            # Eliminiere alle anderen Einträge in dieser Spalte
            for row in range(num_rows):
                if row != current_row and augmented_matrix[row, col] == 1:
                    augmented_matrix[row] = augmented_matrix[row] + augmented_matrix[current_row]
            
            current_row += 1
        
        # Extrahiere die Lösung x1
        x1 = GF2.Zeros(k)
        for i in range(min(k, num_rows)):
            if i < num_rows:
                x1[i] = augmented_matrix[i, -1]

        # Schritt 4: Berechne den Klartext x = x1 * S^-1
        S_inv = np.linalg.inv(self.S)
        plaintext = x1 @ S_inv

        return np.array(plaintext, dtype=int)