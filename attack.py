from mcelice import McElice

"""
McElice Attack Implementation: Information Set Decoding Attack

Die Idee besteht darin, zufällig k Spalten der Generatormatrix G' zuwählen
und zu hoffen, dass diese Positionen des Geheimtextes keine Fehler enthalten.


"""
class InformationSetDecoding(McElice):

    def __init__(self):
        super().__init__()

    def attack(self, ciphertext, I: list[int]) -> str:
        """
        @param I: Die Spalten I in G' sollten nicht fehlerbehaftet sein.
        
        """
        # Implement the attack logic here
        pass