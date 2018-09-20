#! /usr/bin/env python

import time, argparse, logging

import pandas

import bashthebug

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",required=True,help="the csv file downloaded from the Zooniverse containing all the classifcations done to date")
    parser.add_argument("--output_stem",default="test",help="the first part of each of the output files")
    parser.add_argument("--timings",action='store_true',default=False,help="print the time taken for each step")
    options = parser.parse_args()

    # parse the output file to work out the output stem
    output_stem=options.input_file.split("bash-the-bug-classifications")[1].split(".csv")[0]

    # open a log file to record images where the wells cannot be identified
    logging.basicConfig(filename="log/btb-classifications-analyse-"+output_stem+".log",level=logging.INFO,format='%(levelname)s: %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

    print("Reading classifications from CSV file...")
    start=time.time()
    current_classifications=bashthebug.BashTheBugClassifications(zooniverse_file=options.input_file)
    if options.timings:
        print("%.1f seconds" % (time.time()-start))

    current_classifications.extract_classifications()

    current_classifications.create_measurements_table()

    current_classifications.create_users_table()

    start=time.time()
    for sampling_time in ['month','week','day']:

        current_classifications.plot_classifications_by_time(sampling=sampling_time,filename='pdf/'+output_stem+'-classifications-'+sampling_time+'.pdf',add_cumulative=True)
        current_classifications.plot_users_by_time(sampling=sampling_time,filename='pdf/'+output_stem+'-users-'+sampling_time+'.pdf',add_cumulative=True)

    current_classifications.plot_user_classification_distribution(filename="pdf/"+output_stem+'-user-distribution.pdf')

    if options.timings:
        print("%.1f seconds" % (time.time()-start))

    logging.info(current_classifications)

    print("Saving compressed CSV file...")
    start=time.time()
    current_classifications.save_csv("dat/bash-the-bug-classifications.csv.bz2",compression=True)
    if options.timings:
        print("%.1f seconds" % (time.time()-start))

    logging.info(current_classifications.users[["classifications","rank"]][:20])

    # print(current_classifications.gini_coefficient)

    # print("Filtering on study")
    # start=time.time()
    # a.filter_study("CRyPTIC1")
    # print("%.1f seconds" % (time.time()-start))

    # print("Calculating the task duration..")
    # start=time.time()
    # a.calculate_task_durations()
    # print("%.1f seconds" % (time.time()-start))
    #
    # print("Creating misc fields..")
    # start=time.time()
    # a.create_misc_fields()
    # print("%.1f seconds" % (time.time()-start))
    #
    # # can print statistics or plot graphs now
    # pandas.set_option('display.max_columns', 100)
    # print(a.df[:3])




    # print(a.df.task_duration.dtype)
    # print(a.df.viewport_width.dtype)
    # print(a.df.reading_day.dtype)
