
# Introduction
File I/O is an important aspect of data analysis. For this project, we produce two functions, parse_line and read_file which extract data from a text entry and aggregate those findings in dictionary. Together, these functions give the user the ability to specify a threshold for numeric variables and to extract a text value, such as clinical diagnoses, from the entry if the numeric variable is below a threshold. This is a simple way of finding general trends within the data.

This assignment gave us the opportunity to practice working with an unfamiliar file type and the string parsing necessary to extract information from it. Of course, string parsing will be an important skill for bioinformatics in general. We relied on the .split() method, though regex and an even more basic loop approach could have been used.

## Some talking points:

- We did not program defensively and include exception handling for file reading or data types. It is assumed that the VCF file used is not malformed in any way. This might not be a safe assumption and would need to be addressed in future versions of the program.

- We considered several approaches to parsing the string before settling on the split method. The split method, while more complex than regular expressions (regex), is a lower-level approach that is generalizable to many different situations.

- For this assignment, we didn't focus on the meaning of the variables at all. We were very much focused on the code. It is interesting to see how others investigated the meaning of AF_EXAC in the context of clinical diagnosis.

## Files:

project01.py is the main script for this assignment

regex_script.py is a secondary script demonstrating the regex approach to the problem along with different coding conventions

test_script.py uses unittest to assess the correctness of both scripts against different cases and one another

# Pseudocode
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
        if line is legit
            CLNDN values = parse_line(string, threshold)
            add to unique or create key in unique
            
```

# Successes
- Extensive work with .split method and thinking through how to parse strings
- Experience with assigning variables in the constructor and understanding what that does
- Experience with if/then heirarchy and how logic affects code performance

# Struggles
- Git branching and file location was confusing at times
- Incorrect output with some drafts of the code for confusing reasons
- Debugging

# Personal Reflections
## Group Leader (Eric Arnold)
Coding is like playing a musical instrument. It is a highly specialized activity and the user interface is unintuitive. It takes a while to get the hang of, but one you know your scales and can play in tune, you can really jam. This is the struggle before the jam.

## Other member
Other members' reflections on the project

(make sure to mention: positioning the return statement -- I had to change this because it was bugged in the latest version -- and the use of continue as an artifact of dealing with none types)

# Generative AI Appendix
- Used to generate a basic unittest test scheme (It's been a year since I've used it). Model: Sonnet-4.5. Prompt: could you write a basic unittest scheme
- Some syntax checking for regex. I didn't know you could use an f-string directly with regex. Model: Sonnet-4.5. Prompt: how might I write a regex pattern that searches for {keyword}= ; and returns the stuff between the equals and the semicolon


