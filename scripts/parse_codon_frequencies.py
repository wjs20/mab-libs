import re

pat = re.compile(r'([ATGC]{3}\s+\d+\.\d)')

def get_codon_frequencies(raw):
    res = [
        [o.split() for o in pat.findall(line)] 
        for line in raw.replace('U', 'T').splitlines() if line
    ]
    res = dict(list(it.chain.from_iterable(res)))
    res = {codon: float(pct) for codon, pct in res.items()}
    return dict(sorted(res.items())) 

if __name__ == '__main__':
    # https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10029
    raw = """UUU 19.6(  3005)  UCU 16.0(  2450)  UAU 13.1(  2017)  UGU  9.1(  1397)
    UUC 22.0(  3381)  UCC 16.5(  2529)  UAC 16.4(  2519)  UGC 10.3(  1589)
    UUA  6.4(   978)  UCA 10.3(  1577)  UAA  0.6(    93)  UGA  1.2(   177)
    UUG 14.1(  2169)  UCG  3.4(   529)  UAG  0.5(    84)  UGG 13.1(  2012)

    CUU 13.2(  2023)  CCU 16.7(  2563)  CAU 10.2(  1563)  CGU  5.6(   863)
    CUC 18.4(  2818)  CCC 17.0(  2608)  CAC 12.9(  1980)  CGC  9.3(  1429)
    CUA  7.6(  1174)  CCA 15.6(  2388)  CAA 10.3(  1587)  CGA  7.2(  1102)
    CUG 38.8(  5955)  CCG  4.3(   657)  CAG 33.4(  5122)  CGG 10.1(  1558)

    AUU 17.4(  2673)  ACU 14.1(  2172)  AAU 17.4(  2671)  AGU 11.4(  1756)
    AUC 24.8(  3808)  ACC 20.3(  3118)  AAC 21.2(  3248)  AGC 16.4(  2521)
    AUA  6.9(  1053)  ACA 15.7(  2418)  AAA 24.6(  3782)  AGA 10.1(  1557)
    AUG 23.0(  3538)  ACG  4.5(   685)  AAG 38.4(  5895)  AGG 10.2(  1570)

    GUU 11.6(  1780)  GCU 22.4(  3432)  GAU 24.6(  3781)  GGU 12.8(  1968)
    GUC 15.7(  2408)  GCC 25.9(  3973)  GAC 28.1(  4310)  GGC 21.3(  3268)
    GUA  7.8(  1202)  GCA 16.3(  2497)  GAA 28.4(  4355)  GGA 15.8(  2425)
    GUG 30.1(  4628)  GCG  5.0(   765)  GAG 41.1(  6311)  GGG 13.4(  2063)
    """

    f = get_codon_frequencies(raw)
