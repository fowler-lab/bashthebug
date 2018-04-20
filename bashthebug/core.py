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

from tqdm import tqdm

class ZooniverseClassifications():

    def __init__(self,zooniverse_file=None,pickle_file=None):

        self.drug_list={'BDQ':8,'KAN':5,'ETH':6,'AMI':6,'EMB':8,'INH':7,'LEV':7,'MXF':7,'DLM':7,'LZD':7,'CFZ':7,'RIF':7,'RFB':6,'PAS':6}

        if zooniverse_file:

            # read in the raw classifications, parsing the JSON as you go
            self.classifications = pandas.read_csv(zooniverse_file,parse_dates=['created_at'],infer_datetime_format=True,converters={'subject_data':self._parse_json,'metadata':self._parse_json,'annotations':self._parse_json},index_col='classification_id')

            # drop a few of the (unused) columns
            self.classifications.drop(['gold_standard','expert'], axis=1, inplace=True)

            # then create a Boolean column defining if the classification was made during live or not
            self.classifications['live_project']  = [self._get_live_project(q) for q in self.classifications.metadata]

            # create a new dataset containing only live classifications
            self.classifications=self.classifications.loc[self.classifications["live_project"]==True]

            self.classifications['filename']=self.classifications.apply(self._extract_filename2,axis=1)

            self.classifications['plate_image'], self.classifications['drug']=self.classifications['filename'].str.split('-zooniverse-', 1).str

            self.classifications['bashthebug_dilution']=self.classifications.apply(self._parse_annotation,axis=1).astype(int)

            self.classifications["study_id"]=self.classifications.apply(self.determine_study,axis=1)

            self.create_other_tables()

        elif pickle_file:

            # find out the file extension so we can load in the dataset using the right method
            stem, file_extension = os.path.splitext(pickle_file)

            # doing it this way means you can provide either pickle file and it will still work
            self.classifications=pandas.read_pickle(stem+"-classifications"+file_extension)

            self.create_other_tables()

            # self.users=pandas.read_pickle(stem+"-users"+file_extension)
            # self.measurements=pandas.read_pickle(stem+"-measurements"+file_extension)


    def create_other_tables(self):

        # self.classifications.drop(['metadata','annotations','subject_data','filename'], axis=1, inplace=True)

        # create a table of measurements, additional measurements (e.g. Vizion or AMyGDA) can be merged in later
        self.measurements=self.classifications[['plate_image','drug','bashthebug_dilution']].groupby(['plate_image','drug']).agg({'bashthebug_dilution':['median','mean','std','min','max','count']})

        # rename the top level of the columns
        self.measurements.columns.set_levels(['bashthebug_dilution'],level=0,inplace=True)

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

        # self.classifications[['study_id','filename','plate_image','strain','site','duplicate','reader','reading_day','drug']]=self.classifications.apply(self._extract_filename,axis=1)
        self.classifications['filename']=self.classifications.apply(self._extract_filename2,axis=1)

        self.classifications['plate_image'], self.classifications['drug']=self.classifications['filename'].str.split('-zooniverse-', 1).str

        self.classifications['bashthebug_dilution']=self.classifications.apply(self._parse_annotation,axis=1)
        self.classifications['bashthebug_dilution'].astype(int)

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

    def _parse_json(self,data):
        return ujson.loads(data)

    def _get_live_project(self,row):
        try:
            return row['live_project']
        except:
            # apparently some subject metadata doesn't have this? dunno?
            return False

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

    def plot_classifications_by_time(self,sampling=None,colour='#dc2d4c',filename=None,pre_launch=True,add_cumulative=False):

        # create a dataseries of just the date/times
        if not pre_launch:
            a=self.classifications.loc[self.classifications.created_at>datetime.date(2017,4,7)]
            tmp=a[['created_at']]
        else:
            tmp=self.classifications[['created_at']]

        # set it as the index and then re-sample
        tmp.set_index(tmp.created_at,inplace=True)

        self._plot_time_bar(tmp,sampling,colour,filename,pre_launch,add_cumulative,yaxis="Classifications")

    def plot_users_by_time(self,sampling='week',colour='#9ab51e',filename=None,pre_launch=True,add_cumulative=False):

        # create a dataseries of just the date/times
        if not pre_launch:
            a=self.classifications.loc[self.classifications.created_at>datetime.date(2017,4,7)]
            tmp=a[['user_name','created_at']]
        else:
            tmp=self.classifications[['user_name','created_at']]

        foo=tmp[['user_name','created_at']].groupby('user_name').min()
        foo.set_index(foo.created_at,inplace=True)
        tmp=foo.created_at

        self._plot_time_bar(foo,sampling,colour,filename,pre_launch,add_cumulative,yaxis="Users")

    def plot_user_classification_distribution(self,colour="#9ab51e",filename=None):

        stem, file_extension = os.path.splitext(filename)

        # use a square figure
        fig = plt.figure(figsize=(5, 5))
        axes1 = plt.gca()

        axes1.plot(self.users.proportion_user_base,self.users.proportion_total_classifications,color=colour)
        axes1.plot([0,1],[0,1],color=colour,linestyle='dashed')
        axes1.text(0.15,0.65,"Gini-coefficient = %.2f" % self.gini_coefficient,color=colour)
        fig.savefig(stem+"-linear"+file_extension,transparent=True)

        fig = plt.figure(figsize=(5, 5))
        axes1 = plt.gca()

        axes1.set_xscale("log")
        axes1.text(0.0005,0.65,"Gini-coefficient = %.2f" % self.gini_coefficient,color=colour)
        axes1.plot(self.users.proportion_user_base,self.users.proportion_total_classifications,color=colour)
        fig.savefig(stem+"-log"+file_extension,transparent=True)

    def _plot_time_bar(self,data,sampling='week',colour='#e41a1c',filename=None,pre_launch=True,add_cumulative=False,yaxis=None):

        assert sampling in ['week','month','day'], "sampling must be either week, month or day"

        assert filename is not None, 'need to specify a filename with a valid extension'

        if sampling=='week':
            resampled_data=data.resample('W').count()
            bar_width=7.1
        elif sampling=="month":
            resampled_data=data.resample('M').count()
            bar_width=20
        elif sampling=="day":
            resampled_data=data.resample('D').count()
            bar_width=1.01

        resampled_data.columns=['number']

        resampled_data['total']=resampled_data.cumsum()

        fig = plt.figure(figsize=(9, 5))
        axes1 = plt.gca()

        axes1.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))

        # axes1.spines['top'].set_visible(False)
        # axes1.spines['right'].set_visible(False)

        axes1.set_ylabel(yaxis+" per "+sampling,color=colour)
        axes1.xaxis.set_major_locator(mdates.MonthLocator())
        axes1.xaxis.set_major_formatter(mdates.DateFormatter('%b %y'))
        axes1.tick_params('y', colors=colour)
        axes1.bar(resampled_data.index,resampled_data.number,width=bar_width,align='center',lw=0,fc=colour,zorder=10)

        if add_cumulative:
            axes2 = axes1.twinx()
            axes2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
            axes2.tick_params('y', colors='black')
            axes2.xaxis.set_major_locator(mdates.MonthLocator())
            axes2.xaxis.set_major_formatter(mdates.DateFormatter('%b %y'))
            axes2.plot(resampled_data.index,resampled_data.total,zorder=20,color='black')

        fig.savefig(filename,transparent=True)
