"""
Simplified and enhanced interface to optparse, as well as validation methods.
"""

MACROS = dict(
    chrm={ "A" : range(1,20), "S" : ("X","Y"), "F": ("A","X"), "N": ("A","S"), "*" : ("N","M") },
    chrm_plink={ "A" : range(1,20), "X" : (23,), "Y" : (24,), "M" : (26,), "S" : ("X", "Y"), "F": ("A", "X"), "N": ("A","S"), "*" : ("N","M") },
    strain={
        "cccommon": ("129S1", "A_J", "CAST", "NOD", "NZO", "PWK", "WSB"),
        "cc": ("cccommon", "B6"),
        "common": ("cccommon","AKR","BALB","C3H","C57BL","CBA","DBA","LP_J"),
        "sanger": ("common","129P2","129S5","SPRET")
    }
)
