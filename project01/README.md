
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

I don't have much to say other than that because I am a data science student and have to deal with file I/O regularly. It was fun talking through the problem and the challenges associated with python fundamentals.

## Other member

### Connor Crawford

One interesting issue we ran into was dealing with NoneTypes being returned from our parse_line function. Originally, we had the return of an empty list for unsatisfactory conditions positioned after an else statement at the end of the function. This would only return the empty list if the AF_EXAC value wasn't in the line, however, there were situations where we wanted an empty list returned even if AF_EXAC was found in the line. The simple fix ended up being to take away the else statement and just return an empty list, this returned an empty list for any condition other than the ones that met the specific requirements to return a list of diagnosis values.

A language based problem that I ran into came when accessing the values contained in a tuple. It's been a while since I've done much coding in Python, but I have been working in R a lot, where indexing elements starts at 1, not 0. As you can guess Python didn't like when I tried to access the second element of a tuple that only had two values with 2 instead of 1. This assignment definitely helped jog my memory, specifically in regards to string operations, indexing, and working with the control statements in Python.

It's worth noting that there are a few nested for loops in the parse_line function and it does quite a few things. Converting the tuple list into a dictionary to avoid the nested for loops and splitting the parse_line function into smaller functions could make the same "method" of execution more efficient and elegant.

Since it's been quite a while since I worked with Python, this project helped me brush up the basics and get me back on track. Also, I used to believe analysing the genomic data is the hardest part of all, but now I realised that interpreting the data is equally crucial and hard.

### Aaronie Jersha Jenyfred

Since it's been quite a while since I worked with Python, this project helped me brush up the basics and get me back on track. Also, I used to believe analysing the genomic data is the hardest part of all, but now I realised that interpreting the data is equally crucial and hard.

# Generative AI Appendix
- Used to generate a basic unittest test scheme (It's been a year since I've used it). Model: Sonnet-4.5. Prompt: could you write a basic unittest scheme
- Some syntax checking for regex. I don't have extensive experience with regex, and I didn't know you could use an f-string directly with regex. Model: Sonnet-4.5. Prompt: how might I write a regex pattern that searches for {keyword}= ; and returns the stuff between the equals and the semicolon


