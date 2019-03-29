import re


SNP_RE = re.compile(r'(?P<old_allele>[ATCG])>(?P<new_allele>[ATCG])')
SYN_SNP_RE = re.compile(r'(?P<new_allele>[ATCG])=')
INDEL_RE = re.compile(r'.*(?P<new_allele>del.*|ins.*|dup.*)')
PROT_RE = re.compile(r'(?P<aa1>[A-Za-z]{3})(?P<pos>\d+)(?P<aa2>[A-Za-z]{3}|=|\*)')
