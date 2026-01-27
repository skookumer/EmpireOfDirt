from pathlib import Path
import re
from collections import defaultdict

home = Path(__file__).parent
path = home / "clinvar_20190923_short.vcf"


def parse_function(text, threshold_key, retrieval_key, threshold=0.0001):
    '''
    function parses the line and looks for specific keywords using regex, 
    then returns what is between the equal sign and the next semicolon
    '''
    target_pattern = re.compile(rf'{threshold_key}=([^;]+)')                # compiling patterns
    retrieval_pattern = re.compile(rf'{retrieval_key}=([^;]+)')
    if re.search(target_pattern, text):                                     # check if target exists
        value = target_pattern.search(text)
        value = float(value.group(1))
        if value < threshold and re.search(retrieval_pattern, text):        # check if retrieval exists and if value < threshold
            result = retrieval_pattern.search(text)
            result = result.group(1)
            return result.split("|")                                        # return the list split by vertical bars
    return []


def read_file(path, avoid={"not_provided", "not_specified"}):
    '''
    scroll thru the lines of the file and read them one by one
    add unique items to dict and count them
    '''
    unique = defaultdict(int)                                               # use defaultdict to initialize values
    for line in open(path, "r", encoding="utf-8"):
        dn_vals = parse_function(line, "AF_EXAC", "CLNDN")
        for dn in dn_vals:
            unique[dn] += 1
    [unique.pop(key) for key in avoid]                                      # remove unwanted values
    return dict(unique)                                                     # cast as dict for return
    

results = read_file(path)
print(results)