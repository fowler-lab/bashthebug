#! /usr/bin/env python

import dateutil.parser
import datetime
import pandas
import ujson
import numpy
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import bashthebug

class BashTheBugClassifications(bashthebug.ZooniverseClassifications):


    def create_measurements_table(self):

        # self.classifications.drop(['metadata','annotations','subject_data','filename'], axis=1, inplace=True)

        # create a table of measurements, additional measurements (e.g. Vizion or AMyGDA) can be merged in later
        self.measurements=self.classifications[['plate_image','drug','bashthebug_dilution']].groupby(['plate_image','drug']).agg({'bashthebug_dilution':['median','mean','std','min','max','count']})

        # rename the top level of the columns
        self.measurements.columns.set_levels(['bashthebug_dilution'],level=0,inplace=True)

    def create_users_table(self):

        # create a table of users and how many classifications they have done
        self.users=self.classifications[['user_name','created_at']].groupby('user_name').count()

        # rename the total column
        self.users.columns=['classifications']

        # how many users have contributed?
        self.total_users=len(self.users)

        # sort the table so the top users are first
        self.users.sort_values(['classifications'],ascending=False,inplace=True)

        # label each user as whether it is anonymous or not
        self.users['anonymous']=self.users.apply(self._anonymous_user,axis=1)

        # # create a cumulating column
        self.users['cumulative_classifications']=self.users.classifications.cumsum()

        # how many classifications have been done?
        self.total_classifications=self.users.classifications.sum()

        # use this to create a percentage of the total column
        self.users['proportion_total_classifications']=self.users['cumulative_classifications']/self.total_classifications

        # number the users
        self.users['rank']=range(1,self.total_users+1,1)

        # calculate the proportion of the user base
        self.users['proportion_user_base']=self.users['rank']/self.total_users

        # now calculate the Gini coefficient
        area_under_curve=(self.users['proportion_total_classifications'].sum())/self.total_users
        self.gini_coefficient=1-(2*area_under_curve)

    def merge_other_dataset(self,filename=None,new_column=None):

        # find out the file extension so we can load in the dataset using the right method
        stem, file_extension = os.path.splitext(filename)

        assert file_extension in ['.csv','.pkl'], "Only .csv and .pkl file extensions are recognised"

        # read in the datafile in the appropriate way
        if file_extension==".csv":
            other_dataset=pandas.read_csv(filename)
        elif file_extension==".pkl":
            other_dataset=pandas.read_pickle(filename)

        # check that only two columns have been specified
        assert other_dataset.shape[1]==2, "new dataset has more than two columns!"

        # check that the new column is part of thew dataset
        assert new_column in other_dataset.keys(), "specified column "+new_column+" is not in the new dataset!"
        assert 'filename' in other_dataset.keys(), "new dataset does not contain a column named filename to merge on!"

        # split the filename into its components
        other_dataset['plate_image'], other_dataset['drug'] = other_dataset['filename'].str.split('-zooniverse-', 1).str

        # drop the filename column
        other_dataset.drop(['filename'], axis=1, inplace=True)

        # finally create a dataset with a hierarchical index
        summarised_data=other_dataset[['plate_image','drug',new_column]].groupby(['plate_image','drug']).agg({new_column:['median','count']})

        # check that the existing dataframe does not already contain the new column
        assert new_column not in self.measurements.columns.get_level_values(0), "specified column "+new_column+" alreadyx exists in the dataset!"

        # finally, if all that is true, perform the merge
        self.measurements=pandas.merge(self.measurements,summarised_data,left_index=True,right_index=True,how="left")

    def extract_cryptic1_fields(self):

        self.classifications['reading_day']=self.classifications['plate_image'].str.split('-').str[-1].astype(int)
        self.classifications['reader']=self.classifications['plate_image'].str.split('-').str[-2].astype(int)
        self.classifications['replicate']=self.classifications['plate_image'].str.split('-').str[-3].astype(int)
        self.classifications['site']=self.classifications['plate_image'].str.split('-').str[-4].astype(int)

    def determine_study(self,row):

        try:
            if row.filename[:3] in ('H37','CRY'):
                return "CRyPTIC1"
            else:
                return "CRyPTIC2"
        except:
            return "Unknown"

    def extract_classifications(self):

        self.drug_list={'BDQ':8,'KAN':5,'ETH':6,'AMI':6,'EMB':8,'INH':7,'LEV':7,'MXF':7,'DLM':7,'LZD':7,'CFZ':7,'RIF':7,'RFB':6,'PAS':6}

        self.classifications['filename']=self.classifications.apply(self._extract_filename2,axis=1)

        self.classifications['plate_image'], self.classifications['drug']=self.classifications['filename'].str.split('-zooniverse-', 1).str

        self.classifications['bashthebug_dilution']=self.classifications.apply(self._parse_annotation,axis=1).astype(int)

        self.classifications["study_id"]=self.classifications.apply(self.determine_study,axis=1)



    def calculate_consensus_median(self):

        # create a consensus based on the median
        self.consensus_median=self.classifications[["filename","bashthebug_dilution"]].groupby('filename').median()

        # rename it
        self.consensus_median.columns=['bashthebug_median']

        # merge it back into the dataset
        self.classifications=pandas.merge(self.classifications,self.consensus_median,left_on="filename",how='left',right_index=True)

        # calculate for each classification how far it is away from the consensus
        self.classifications['median_delta']=self.classifications['bashthebug_dilution']-self.classifications['bashthebug_median']

    def calculate_task_durations(self):

        self.classifications['task_duration']=self.classifications.apply(self._task_duration,axis=1)

    def filter_study(self,study):

        self.classifications=self.classifications.loc[self.classifications["study_id"]==study]

    def filter_readingday(self,reading_day):

        self.classifications=self.classifications.loc[self.classifications["reading_day"]==reading_day]

    def create_misc_fields(self):

        self.classifications['user_language']=self.classifications.apply(self._user_language,axis=1)
        self.classifications['viewport_width']=self.classifications.apply(self._parse_viewport_width,axis=1)
        self.classifications['viewport_height']=self.classifications.apply(self._parse_viewport_height,axis=1)

    def save_pickle(self,filename):

        stem, file_extension = os.path.splitext(filename)
        self.classifications.to_pickle(stem+"-classifications"+file_extension)
        self.users.to_pickle(stem+"-users"+file_extension)
        self.measurements.to_pickle(stem+"-measurements"+file_extension)

    def save_csv(self,filename):

        stem, file_extension = os.path.splitext(filename)
        self.classifications.to_csv(stem+"-classifications"+file_extension)
        self.users.to_csv(stem+"-users"+file_extension)
        self.measurements.to_csv(stem+"-measurements"+file_extension)

    def __repr__(self):

        line="%30s %7i\n" % ("Total classifications:",self.total_classifications)
        line+="%30s %7i\n" % ("Total users:",self.total_users)
        line+="%30s %7.2f\n" % ("Gini coefficient:",self.gini_coefficient)
        line+="\n"
        top_10=(100*self.users.classifications[:10].sum())/self.total_classifications
        line+="%30s %7.1f %%\n" % ("Top   10 users have done:",top_10)
        top_100=(100*self.users.classifications[:100].sum())/self.total_classifications
        line+="%30s %7.1f %%\n" % ("Top  100 users have done:",top_100)
        top_1000=(100*self.users.classifications[:1000].sum())/self.total_classifications
        line+="%30s %7.1f %%\n" % ("Top 1000 users have done:",top_1000)

        return(line)


    def _extract_filename2(self,row):
        try:
            for i in row.subject_data[str(row.subject_ids)]:
                if (".png" in i) or i=="Filename":
                    return(row.subject_data[str(row.subject_ids)][i][:-4])
        except:
            print("Problem parsing "+row.classification_id)


    def _extract_filename(self,row):
        strain=None
        site=None
        duplicate=None
        reader=None
        reading_day=None
        drug=None
        study_id=None
        filename=None
        plate_image=None
        try:
            for i in row.subject_data[str(row.subject_ids)]:
                if (".png" in i) or i=="Filename":
                    filename=row.subject_data[str(row.subject_ids)][i][:-4]
                    try:
                        if filename[:3] in ('H37','CRY'):
                            study_id="CRyPTIC1"
                        else:
                            study_id="Unknown"
                    except:
                        study_id="Unknown"
                    tmp=filename.split('-zooniverse-')
                    plate_image=tmp[0]
                    drug=tmp[1].split('.')[0]
                    foo=tmp[0].split('-')
                    if foo[0]=='CRY':
                        strain=foo[0]+"-"+foo[1]
                        site=foo[2]
                        duplicate=foo[3]
                        reader=foo[4]
                        reading_day=foo[5]
                    elif foo[0]=="H37rV":
                        strain=foo[0]
                        site=foo[1]
                        duplicate=foo[2]
                        reader=foo[3]
                        reading_day=foo[4]

        except:
            print("Problem parsing "+row.classification_id)
        return(pandas.Series([study_id,filename,plate_image,strain,site,duplicate,reader,reading_day,drug]))

    def _task_duration(self,row):
        try:
            start=(dateutil.parser.parse(row.metadata['started_at']))
            end=(dateutil.parser.parse(row.metadata['finished_at']))
            duration = (end-start)/pandas.Timedelta(1, unit='s')
            return(duration)
        except:
            return(pandas.Timedelta(0, unit='s'))

    def _user_language(self,row):
        try:
            return row.metadata['user_language']
        except:
            return ""

    def _parse_viewport_width(self,row):
        return row.metadata['viewport']['height']

    def _parse_viewport_height(self,row):
        return row.metadata['viewport']['height']

    def _anonymous_user(self,row):
        if row.name[0:13]=="not-logged-in":
            return True
        else:
            return False

    def _parse_annotation(self,row):
        try:
            answer_text=row.annotations[0]["value"]
            if ("No Growth in either" in answer_text) or ("No Growth in one" in answer_text):
                return -1
            elif ("No Growth in wells" in answer_text) or ("No Growth in all" in answer_text):
                return 1
            elif ("Growth in all" in answer_text):
                return int(self.drug_list[row.drug]+1)
            elif "Cannot classify" in answer_text:
                return 0
            elif ("dose" in answer_text) or ("identify" in answer_text):
                try:
                    return int(row.annotations[1]["value"])
                except:
                    return 0
            else:
                return 0
        except:
            return 0
