#! /usr/bin/env python

import time
import argparse

import pandas

import bashthebug

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",required=True,help="the csv file downloaded from the Zooniverse containing all the classifcations done to date")
    parser.add_argument("--output_stem",default="test",help="the first part of each of the output files")
    parser.add_argument("--timings",action='store_true',default=False,help="print the time taken for each step")
    options = parser.parse_args()

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

        current_classifications.plot_classifications_by_time(sampling=sampling_time,filename='pdf/'+options.output_stem+'-classifications-'+sampling_time+'.pdf',pre_launch=True,add_cumulative=True)
        current_classifications.plot_users_by_time(sampling=sampling_time,filename='pdf/'+options.output_stem+'-users-'+sampling_time+'.pdf',pre_launch=True,add_cumulative=True)

    current_classifications.plot_user_classification_distribution(filename="pdf/"+options.output_stem+'-user-distribution.pdf')

    if options.timings:
        print("%.1f seconds" % (time.time()-start))

    print(current_classifications)

    print("Saving PKL file...")
    start=time.time()
    current_classifications.save_pickle("dat/bash-the-bug-classifications-"+options.output_stem+".pkl")
    if options.timings:
        print("%.1f seconds" % (time.time()-start))


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
