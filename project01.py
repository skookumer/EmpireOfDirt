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
                        clndn = tup2[1].split(",")

                        # now we can define a diagnosis list to store all the diagnoses of interest
                        diagnosis_list = []
                        # loop over each diagnosis
                        for diagnosis in clndn:
                            # make sure that the diagnosis isn't one we don't want
                            if diagnosis != 'not_specified' and diagnosis != 'not_provided':
                                # add it to the diagnosis list
                                diagnosis_list.append(diagnosis)

                        return diagnosis_list

            # if the AF_EXAC value doesn't meet the required threshold, return an empty list
            else:
                return []
