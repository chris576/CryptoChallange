"""
Bitmap Implementation Package für McElice Kryptosystem
Optimiert mit gmpy2 für schnelle Bit-Operationen
"""

from .mcelice_bitmap import (
    BitMatrix, 
    BitVector, 
    McElice, 
    parity_check_matrix, 
    decode_syndrome_bitmap,
    popcount,
    bit_and,
    bit_or,
    bit_xor,
    bit_not,
    test_bit,
    set_bit,
    clear_bit
)
from .isd_attack import InformationSetDecoding

__all__ = [
    'BitMatrix',
    'BitVector', 
    'McElice',
    'parity_check_matrix',
    'decode_syndrome_bitmap',
    'InformationSetDecoding',
    'popcount',
    'bit_and',
    'bit_or',
    'bit_xor',
    'bit_not',
    'test_bit',
    'set_bit',
    'clear_bit'
]
