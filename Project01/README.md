
# Introduction
File I/O is an important aspect of data analysis. For this project, we produce two functions, parse_line and read_file which extract data from a text entry and aggregate those findings in dictionary. Together, these functions give the user the ability to specify a threshold for numeric variables and to extract a text value, such as clinical diagnoses, from the entry if the numeric variable is below a threshold. This is a simple way of finding general trends within the data.

This assignment gave us the opportunity to practice working with an unfamiliar file type and the string parsing necessary to extract information from it. Of course, string parsing will be an important skill for bioinformatics in general. We relied on the .split() method, though regex and an even more basic loop approach could have been used.

# Pseudocode
Put pseudocode in this box:

```python

function parse_line(string, threshold):
    split_string = string split on semicolon
    if "AF_EXEC" in split_string:
        get AF_EXEC value
        if value < thresold:
            get CLNDN values
            return CLNDN values

function read_file(path, threshold):
    unique is a dict
    for line in open(path)
        if line is legitimate
            CLNDN values = parse_line(string, threshold)
            add to unique or create key in unique
            
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
