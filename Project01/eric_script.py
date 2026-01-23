from pathlib import Path
import re

home = Path(__file__).parent
path = home / "clinvar_20190923_short.vcf"


def parse_function(text, threshold_key, retrieval_key, threshold=0.0001):
    target_pattern = re.compile(rf'{threshold_key}=([^;]+)')
    retrieval_pattern = re.compile(rf'{retrieval_key}=([^;]+)')
    if re.search(target_pattern, text):
        value = target_pattern.search(text)
        value = float(value.group(1))
        if value < threshold and re.search(retrieval_pattern, text):
            result = retrieval_pattern.search(text)
            result = result.group(1)
            return result.split("|")
    return []


def read_file(path, avoid = {"not_provided", "not_specified"}):
    unique = {}
    for line in open(path, "r", encoding="utf-8"):
        dn_vals = parse_function(line, "AF_EXAC", "CLNDN")
        if len(dn_vals) > 0:
            for dn in dn_vals:
                if dn not in avoid:
                    if dn in unique:
                        unique[dn] += 1
                    else:
                        unique[dn] = 1
    return unique
    

results = read_file(path)
print(results)