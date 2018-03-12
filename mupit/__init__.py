
# define the VEP consequence types for loss of function and missense variants
LOF_CQ = set(["stop_gained", "splice_acceptor_variant", "splice_donor_variant",
    "frameshift_variant", "initiator_codon_variant", "start_lost",
    "conserved_exon_terminus_variant"])
MISSENSE_CQ = set(["missense_variant", "stop_lost", "inframe_deletion",
    "inframe_insertion", "coding_sequence_variant",
    "protein_altering_variant"])
SYNONYMOUS_CQ = set(["synonymous_variant"])
