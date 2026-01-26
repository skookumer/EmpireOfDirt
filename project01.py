def parse_line(clinvar_line, threshold):

    # define the accepted data (allowed_data = 'AF_EXAC')
    allowed_data = 'AF_EXAC'

    # determine if AF_EXAC is present
    if allowed_data in clinvar_line:
        # if AF_EXAC is present, split the line by semicolons to separate all information variables in the data
        semicol_split = clinvar_line.split(";")

        # create a tuple list, this will be used to add tuples that package info variable descriptions and their values
        tuple_list = []
        # loop through the range of values that make up the length of the semicol_split list
        for i in range(len(semicol_split)):
            # split each element of the semicol_split list by "=" and pack the two values into a tuple
            equals_split = tuple(semicol_split[i].split("="))
            # append the tuple to the tuple list
            tuple_list.append(equals_split)

        # now that we have a tuple for every info variable, and it's value, we can loop through the list and apply our thresholds and parameters
        for tup in tuple_list:
            # find the tuple with the AF_EXAC info and check that the associated value meets the required threshold
            if allowed_data in tup[0] and float(tup[1]) < threshold:

                # now that we've verified the rarity value meets the required threshold we need to extract out the associated clinical diagnosis with that variant
                # first we define the info variable we want to identify - 'CLNDN'
                clin_diagnosis = 'CLNDN'
                # now we can loop through the tuple list to identify which tuple has the info we need
                for tup2 in tuple_list:
                    # check to see if we've found the clinical diagnosis tuple
                    if clin_diagnosis in tup2[0]:
                        # if we have, split associated diagnosis of the tuple by commas to create a list of associated diagnoses
                        clndn = tup2[1].split("|")

                        # now we can define a diagnosis list to store all the diagnoses of interest
                        diagnosis_list = []
                        # loop over each diagnosis
                        for diagnosis in clndn:
                            # make sure that the diagnosis isn't one we don't want
                            if 'not_specified' not in diagnosis or 'not_provided' not in diagnosis:
                                # add it to the diagnosis list
                                diagnosis_list.append(diagnosis)

                        return diagnosis_list

            # if the AF_EXAC value doesn't meet the required threshold, return an empty list
            else:
                return []
                
def read_file(vcf_file, threshold = 0.0001):
    # create a dictionary to count the number of times a specific disease is observed
    disease_count = {}

    # opens the file
    with open(vcf_file, 'r') as file:

        # reads the file line by line
        for line in file:

            # skips the description part of the vcf file
            if line.startswith('#'):
                continue

            # to pass the line to parse line function to get the disease list
            disease_list = parse_line(line, threshold)

            # to skip the line if None or an empty list is returned
            if not disease_list:
                continue

            # create a loop for the returned list
            for disease in disease_list:
                # to check if the disease is present in the dictionary
                if disease in disease_count:
                    # if present, increment the count by 1
                    disease_count[disease] += 1
                else:
                    # if not present, add the disease to the dictionary with a count of 1
                    disease_count[disease] = 1

    return disease_count

if __name__ == "__main__":
    pprint(read_file("clinvar_20190923_short.vcf"))
