from mcelice import McElice
import numpy as np
import galois

"""
McElice Attack Implementation: Information Set Decoding Attack

Die Idee besteht darin, zufällig k Spalten der Generatormatrix G' zuwählen
und zu hoffen, dass diese Positionen des Geheimtextes keine Fehler enthalten.


"""
class InformationSetDecoding(McElice):

    def __init__(self):
        super().__init__()
        self.GF2 = galois.GF(2)

    def hamming_weight(self, v):
        """Hamming-Gewicht eines Vektors (numpy oder galois)"""
        # galois-Elemente nach numpy ints
        v = np.array(v, dtype=np.int8)
        return int(np.sum(v))
    
    def G_pub(self, G) -> 'InformationSetDecoding':
        self.G_pub = G
        return self
    
    def check_candidate(self, cipher, candidate, G_pub, t):
        """
        Prüft, ob ein rekonstruierter Kandidat korrekt ist.

        y:            empfangenes Wort (Länge n)
        x_candidate:  Kandidat für Klartext (Länge k)
        G_pub:        öffentliche Generatormatrix G' (k x n)
        t:            erwartete Fehleranzahl
        Rückgabe: (bool, int) - ob korrekt, und das Hamming-Gewicht des Fehlers
        """
        y_gf  = self.GF2(cipher)
        x_gf  = self.GF2(candidate)
        G_gf  = self.GF2(G_pub)

        # c' = x' * G'
        c_prime = x_gf @ G_gf

        # r = y - c'  (in GF(2): XOR)
        r = y_gf + c_prime

        # Fehlergewicht prüfen
        wt = self.hamming_weight(r)
        return wt == t, wt        

        
    def attack(self, ciphertext) -> str:
        from lgs import invert
        """
        1. Wähle zufällig k Positionen aus den n Positionen der Generatormatrix G'.
        2. G' mit k*k Spalten als fehlerfrei annehmen und invertieren. Falls G' an diesen Positionen nicht invertierbar ist, wähle neue Positionen.
        3. Entschlüssele den Geheimtext an diesen Positionen.
        4. Berechne den Fehlervektor und überprüfe das Hamming-Gewicht.
        5. Wiederhole den Vorgang, bis der Fehlervektor das erwartete Fehlergewicht hat.
        """
        # Implement the attack logic here
        self.k = self.G_pub.shape[0]
        self.n = self.G_pub.shape[1]
        # Fehlerkapazität: falls nicht gesetzt, Standard 1 (z.B. Hamming)
        t = 0
        max_attempts = 1000000  # Maximale Anzahl an Versuchen
        
        while max_attempts > 0:
            max_attempts -= 1
            
            # 1. Zufällig k Positionen auswählen
            positions = np.random.choice(self.n, self.k, replace=False)
            G_sub = self.G_pub[:, positions]
            
            # 2. Prüfen, ob G_sub invertierbar ist
            try:
                G_inv = invert(G_sub)
            except np.linalg.LinAlgError:
                continue  # Nicht invertierbar, neue Positionen wählen
            
            # 3. Entschlüsseln des Geheimtextes an diesen Positionen
            y_sub = self.GF2(ciphertext)[positions]
            x_candidate = self.GF2(G_inv) @ y_sub
            
            # 4. Fehlervektor berechnen und Hamming-Gewicht prüfen
            is_correct, wt = self.check_candidate(ciphertext, x_candidate, self.G_pub, t)
            print(f"Versuch {1000000 - max_attempts}, Fehlergewicht: {wt}")
            if is_correct:
                return ''.join(map(str, x_candidate.astype(np.int8)))
            
        return "Angriff fehlgeschlagen: Maximale Versuche überschritten."