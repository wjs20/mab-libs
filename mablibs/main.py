import functools

from mablibs import codons, optimization
from mablibs.mutagenesis import Mutagenesis
from mablibs.strategies import RandomizationStrategy
from mablibs.templates import Template


def to_nucleotides(seq_aa):
    return "".join([codons.AA2CODON[aa][0] for aa in seq_aa])


def main():
    # http://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/therasummary/?status=Active&yearprop=2017&format=Whole+mAb+ADC&dev_tech=na&clintrial=Approved&struc95to98=None&cond_disc=na&target=TNFRSF17&isotype=G1&notes=&cond_approved=Multiple+myeloma&companies=GlaxoSmithKline&light2=na&yearrec=2018&INN=belantamab&light1=DIQMTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKLLIYYTSNLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYRKLPWTFGQGTKLEIK&struc99=None&cond_active=na&heavy1=QVQLVQSGAEVKKPGSSVKVSCKASGGTFSNYWMHWVRQAPGQGLEWMGATYRGHSDTYYNQKFKGRVTITADKSTSTAYMELSSLRSEDTAVYYCARGAIYDGYDVLDNWGQGTLVTVSS&heavy2=na&struc100=None
    nums = "1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112A	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128"
    aas = "Q	V	Q	L	V	Q	S	G	A	-	E	V	K	K	P	G	S	S	V	K	V	S	C	K	A	S	G	G	T	F	-	-	-	-	S	N	Y	W	M	H	W	V	R	Q	A	P	G	Q	G	L	E	W	M	G	A	T	Y	R	G	-	-	H	S	D	T	Y	Y	N	Q	K	F	K	-	G	R	V	T	I	T	A	D	K	S	T	S	T	A	Y	M	E	L	S	S	L	R	S	E	D	T	A	V	Y	Y	C	A	R	G	A	I	Y	D	G	Y	D	V	L	D	N	W	G	Q	G	T	L	V	T	V	S	S"
    belantamab_VH_numbered = dict(
        enumerate(zip(nums.split(), aas.split()))
    )  # HCDR3 is at [104:117]
    belantamab_VH_seq = aas.replace("	", "")
    hcdr3_aa = belantamab_VH_seq[100:120]
    hcdr3_nucleotides = to_nucleotides(hcdr3_aa)
    template = Template(hcdr3_nucleotides)
    substitutions = set(template.amino_acids) - {"C", "M"}
    prm = dict(zip(range(0, 16, 2), [list(substitutions)] * 8))
    randomization = RandomizationStrategy(prm, 4)
    restriction_sites = optimization.compile_restriction_site_regex(
        "NheI", "NotI", "XhoI", "NcoI", "DraI"
    )
    nucleotide_repeats = optimization.compile_nucleotide_repeat_regex()
    constraints = [
        functools.partial(
            optimization.pattern_not_found, patterns=restriction_sites
        ),
        functools.partial(
            optimization.is_below_gc_content_threshold, threshold=0.75
        ),
        optimization.is_not_palindromic,
    ]

    species = "e_coli"
    optimizer = optimization.DNAOptimizer(constraints, species)
    mutagenesis = Mutagenesis(
        randomization, template, species, optimizer=optimizer
    )
    for t in mutagenesis.generate_library():
        t.nucleotides


if __name__ == "__main__":
    exit(main())
