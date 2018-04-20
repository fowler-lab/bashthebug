#! /usr/bin/env python

import bashthebug
import time
import pandas

input_file="../download/bash-the-bug-classifications-2018-04-07.csv"
a=input_file.split('.')[-2].split("-")
date=a[-3]+"-"+a[-2]+"-"+a[-1]

print("Instantiating class...")
start=time.time()
current_classifications=bashthebug.ZooniverseClassifications(zooniverse_file=input_file)
print("%.1f seconds" % (time.time()-start))

print("Plotting graphs...")
start=time.time()
for sampling_time in ['month','week','day']:

    current_classifications.plot_classifications_by_time(sampling=sampling_time,filename='pdf/'+date+'-classifications-'+sampling_time+'.pdf',pre_launch=True,add_cumulative=True)
    current_classifications.plot_users_by_time(sampling=sampling_time,filename='pdf/'+date+'-users-'+sampling_time+'.pdf',pre_launch=True,add_cumulative=True)

current_classifications.plot_user_classification_distribution(filename="pdf/"+date+'-user-distribution.pdf')
print("%.1f seconds" % (time.time()-start))

print(current_classifications)

print("Saving pkl...")
start=time.time()
current_classifications.save_pickle("dat/bash-the-bug-classifications-"+date+".pkl")
print("%.1f seconds" % (time.time()-start))

# print("Saving csv...")
# start=time.time()
# current_classifications.save_csv("dat/bash-the-bug-classifications-"+date+".csv")
# print("%.1f seconds" % (time.time()-start))

# print("Parsing filename and extracting classifications..")
# start=time.time()
# a.extract_classifications()
# print("%.1f seconds" % (time.time()-start))

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
