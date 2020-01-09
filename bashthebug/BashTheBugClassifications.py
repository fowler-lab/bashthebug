#! /usr/bin/env python

import os, math

import pandas, numpy
from tqdm import tqdm

import pyniverse

class BashTheBugClassifications(pyniverse.Classifications):

    def _remove_values_from_list(self,the_list,threshold):
        return numpy.array([value for value in the_list if value >= threshold]).astype(int)


    def _custom_aggregate_classifications(self,series):

        n_total=None
        n_failed=None
        n_cannot_read=None
        n_valid=None
        median=None
        mean=None
        std=None
        mmin=None
        mmax=None

        # store the classifications as a list
        classifications=numpy.array(series).astype(int)

        # how many classifications do we have?
        count=len(classifications)

        n_failed=numpy.sum(classifications<-2)

        # now count how many 'cannot read' codes there are
        n_cannot_read = numpy.sum((classifications!=-1) & (classifications!=-2))

        n_valid=numpy.sum(classifications>0)

        # first check we have enough samples
        if len(classifications)>10:

            # now filter out failed classifications
            classifications = self._remove_values_from_list(classifications,-2)

            # if over half the volunteers have said they cannot read the image, return cannot read
            proportion_failed = n_cannot_read/len(classifications)

            if proportion_failed>=0.5 or n_valid<=5:

                failed_code=-1
                median=failed_code
                mean=failed_code
                std=failed_code
                mmin=failed_code
                mmax=failed_code

            else:

                # filter out the cannot read codes
                classifications = self._remove_values_from_list(classifications,1)

                # now finally we can apply the median
                median=math.ceil(numpy.median(classifications))

                mean=numpy.mean(classifications)

                std=numpy.std(classifications)

                mmin=numpy.min(classifications)

                mmax=numpy.max(classifications)

        return(count,n_failed,n_cannot_read,n_valid,median,mean,std,mmin,mmax)

    def create_measurements_table(self,index='PLATEIMAGE'):

        assert index in ['PLATEIMAGE','PLATE'], 'specified index not recognised!'

        # create a table of measurements, additional measurements (e.g. Vizion or AMyGDA) can be merged in later
        if index=='PLATEIMAGE':
            # self.measurements=self.classifications[['plate_image','drug','bashthebug_dilution']].groupby(['plate_image','drug']).agg({'bashthebug_dilution':['median','mean','std','min','max','count']})
            foo=self.classifications[['plate_image','drug','bashthebug_dilution']].groupby(['plate_image','drug']).agg(self._custom_aggregate_classifications)

            self.measurements=pandas.DataFrame(foo['bashthebug_dilution'].tolist(),index=foo.index)

        else:
            # self.measurements=self.classifications[['plate','reading_day','drug','bashthebug_dilution']].groupby(['plate','reading_day','drug']).agg({'bashthebug_dilution':['median','mean','std','min','max','count']})
            foo=self.classifications[['plate','reading_day','drug','bashthebug_dilution']].groupby(['plate','reading_day','drug']).agg(self._custom_aggregate_classifications)

            self.measurements=pandas.DataFrame(foo['bashthebug_dilution'].tolist(),index=foo.index)

        self.measurements.columns=['count','n_failed','n_cannot_read','n_valid','median','mean','std','min','max']
        # self.measurements.columns=['count','median','mean','std','min','max']
        self.measurements=self.measurements[['median','mean','std','min','max','count','n_failed','n_cannot_read','n_valid']]
        # self.measurements=self.measurements[['median','mean','std','min','max','count']]

        # self.classifications.drop(['metadata','annotations','subject_data','filename'], axis=1, inplace=True)

        # rename the top level of the columns
        # self.measurements.columns.set_levels(['bashthebug_dilution'],level=0,inplace=True)

    def create_durations_table(self,index='PLATEIMAGE'):

        assert 'task_duration' in self.classifications.columns, "task_duration not in CLASSIFICATIONS table; make sure you have run the calculate_task_durations() method!"

        assert index in ['PLATEIMAGE','PLATE'], 'specified index not recognised!'

        if index=='PLATEIMAGE':
            self.durations=self.classifications[['plate_image','drug','task_duration']].groupby(['plate_image','drug']).agg({'task_duration':['mean','std']})
        else:
            self.durations=self.classifications[['plate','reading_day','drug','task_duration']].groupby(['plate','reading_day','drug']).agg({'task_duration':['median','mean','std','min','max','count']})

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

        # self.classifications['reading_day']=self.classifications['plate_image'].str.split('-').str[-1].astype(int)
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

    def extract_reading_day(self,row):

        if row['study_id']=='CRyPTIC1':
            reading_day=int(row['plate_image'].split('-')[-1])
        elif row['study_id']=='CRyPTIC2':
            if 'UKMYC' in row['plate_image']:
                reading_day=int(row['plate_image'].split("-")[-2])
            else:
                reading_day=int(row['plate_image'].split("-")[-1])
        else:
            reading_day=None
        return(reading_day)

    def extract_site(self,row):

        if row['study_id']=='CRyPTIC1':
            site=row['plate_image'].split('-')[-4]
        elif row['study_id']=='CRyPTIC2':
            site=row['plate_image'][:2]
        else:
            site=None
        return(site)

    def _extract_plateimage(self,row):

        # find out the filename
        filename=None
        try:
            for i in row.subject_data[str(row.subject_ids)]:
                if (".png" in i) or (".jpg" in i) or i in ["Filename","Image"]:
                    filename=row.subject_data[str(row.subject_ids)][i][:-4]
        except:
            print("Problem parsing "+row.classification_id)

        # find out the study id
        if filename==None:
            study_id=None
        elif filename[:3] in ('H37','CRY'):
            study_id="CRyPTIC1"
        else:
            study_id="CRyPTIC2"

        # extract plate
        plate=None
        if 'plate_image' in row.keys() and row['plate_image'] is not None:
            if "UKMYC" in row['plate_image']:
                foo=row['plate_image'][:-7]
            else:
                foo=row['plate_image']
            location=foo.rfind("-")
            plate=foo[:location]

        # extract drug
        drug=None
        if filename is not None:
            drug=filename[-3:]

        # extract plate_image and plate_design
        if filename is None:
            plate_image=None
            plate_design=None
        elif "UKMYC" in filename:
            plate_image=filename.split("-UKMYC")[0]
            plate_design=filename.split(plate_image+"-")[1].split("-zooniverse")[0]
        else:
            if self.flavour=='regular':
                plate_image=filename.split('-zooniverse-')[0]
            elif self.flavour=='pro':
                plate_image=filename.split('-discrepancy-')[0]
            plate_design="UKMYC5"

        # extract site and reading_day
        if study_id=='CRyPTIC1':
            reading_day=int(plate_image.split('-')[-1])
            site=plate_image.split('-')[-4]
        elif study_id=='CRyPTIC2':
            site=plate_image[:2]
            if 'UKMYC' in plate_image:
                reading_day=int(plate_image.split("-")[-2])
            else:
                reading_day=int(plate_image.split("-")[-1])
        else:
            site=None
            reading_day=None

        return(pandas.Series([filename,plate_image,plate_design,drug,plate,study_id,reading_day,site]))

    def extract_classifications(self,flavour='regular'):

        assert flavour in ['regular','pro'], "unrecognised flavour! "

        # save the flavour
        self.flavour=flavour

        # tqdm.pandas(desc='extracting filename')
        # self.classifications['filename']=self.classifications.progress_apply(self._extract_filename2,axis=1)

        tqdm.pandas(desc='extracting metadata')
        self.classifications[['filename','plate_image','plate_design','drug','plate','study_id','reading_day','site']]=self.classifications.progress_apply(self._extract_plateimage,axis=1)

        # tqdm.pandas(desc='extracting drug')
        # self.classifications['drug']=self.classifications.progress_apply(self._extract_drug,axis=1)
        #
        # tqdm.pandas(desc='extracting plate')
        # self.classifications['plate']=self.classifications.progress_apply(self._extract_plate,axis=1)

        tqdm.pandas(desc='calculating dilution')
        self.classifications['bashthebug_dilution']=self.classifications.progress_apply(self._parse_annotation,axis=1).astype(int)

        # tqdm.pandas(desc='extracting study')
        # self.classifications["study_id"]=self.classifications.progress_apply(self.determine_study,axis=1)

        # tqdm.pandas(desc='extracting reading day')
        # self.classifications['reading_day']=self.classifications.progress_apply(self.extract_reading_day,axis=1)

        # tqdm.pandas(desc='extracting site')
        # self.classifications['site']=self.classifications.progress_apply(self.extract_site,axis=1)

    def calculate_consensus_median(self):

        # create a consensus based on the median
        self.consensus_median=self.classifications[["filename","bashthebug_dilution"]].groupby('filename').median()

        # rename it
        self.consensus_median.columns=['bashthebug_median']

        # merge it back into the dataset
        self.classifications=pandas.merge(self.classifications,self.consensus_median,left_on="filename",how='left',right_index=True)

        # calculate for each classification how far it is away from the consensus
        self.classifications['median_delta']=self.classifications['bashthebug_dilution']-self.classifications['bashthebug_median']

    def filter_study(self,study):

        self.classifications=self.classifications.loc[self.classifications["study_id"]==study]

        self.total_classifications=len(self.classifications)

    def filter_readingday(self,reading_day):

        self.classifications=self.classifications.loc[self.classifications["reading_day"]==reading_day]

        self.total_classifications=len(self.classifications)

    def _extract_filename2(self,row):
        try:
            for i in row.subject_data[str(row.subject_ids)]:
                if (".png" in i) or (".jpg" in i) or i in ["Filename","Image"]:
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


    def _parse_annotation(self,row):

        if 'task_label' not in row.annotations[0]:
            return(-100)

        else:

            # First thing is to work out what the question/task structure is
            task_label=row.annotations[0]["task_label"]

            if "being mindful" in task_label:
                question_type="pro_v1"
            elif "Having looked" in task_label:
                question_type="regular_v1"
            elif "stops, please choose the number" in task_label:
                question_type="regular_v2"
            elif "choose the number" in task_label:
                question_type="regular_v2"
            elif "Mark the first well contain" in task_label:
                question_type="testing"
            else:
                raise ValueError("cannot determine type of task:"+ task_label)

            # ignore the classifications where a red cross is placed
            if question_type=="testing":

                return(-101)

            else:

                answer_text=row.annotations[0]["value"]

                if row['plate_design']=="UKMYC5":
                    drug_list={'BDQ':8,'KAN':5,'ETH':6,'AMI':6,'EMB':8,'INH':7,'LEV':7,'MXF':7,'DLM':7,'LZD':7,'CFZ':7,'RIF':7,'RFB':6,'PAS':6}
                elif row['plate_design']=="UKMYC6":
                    drug_list={'BDQ':8,'KAN':5,'ETH':6,'AMI':7,'EMB':8,'INH':10,'LEV':7,'MXF':7,'DLM':7,'LZD':7,'CFZ':7,'RIF':9,'RFB':6}
                else:
                    raise ValueError("plate design not found")

                if answer_text is None:
                    return(-102)

                elif question_type=='regular_v1':

                    if ("No Growth in either" in answer_text) or ("No Growth in one" in answer_text):
                        return -2
                    elif ("No Growth in wells" in answer_text) or ("No Growth in all" in answer_text):
                        return 1
                    elif ("Growth in all" in answer_text):
                        return int(drug_list[row.drug]+1)
                    elif "Cannot classify" in answer_text:
                        return -1
                    elif len(row.annotations)>1 and row.annotations[1]["value"] is not None:
                        try:
                            return int(row.annotations[1]["value"])
                        except:
                            return(-103)
                    else:
                        return(-104)

                elif question_type=="regular_v2":

                    if ("No Growth in either" in answer_text) or ("No Growth in one" in answer_text):
                        return -2
                    elif ("No Growth in wells" in answer_text) or ("No Growth in all" in answer_text):
                        return 1
                    elif ("Growth in all" in answer_text):
                        return int(drug_list[row.drug]+1)

                    elif "Cannot classify" in answer_text:
                        return -1
                    elif row.annotations[0]["value"].isnumeric():
                        return int(row.annotations[0]["value"])
                    else:
                        return(-105)

                elif question_type=="pro_v1":

                    if ("No Growth in either" in answer_text) or ("No Growth in one" in answer_text):
                        return -2
                    elif ("No Growth in wells" in answer_text) or ("No Growth in all" in answer_text):
                        return 1
                    elif ("Growth in all" in answer_text):
                        return int(self.drug_list[row.drug]+1)
                    elif "Cannot classify" in answer_text:
                        if row.annotations[1]["value"]=="Skip wells":
                            return -10
                        elif row.annotations[1]["value"]=="Trailing pattern":
                            return -11
                        elif row.annotations[1]["value"]=="Contamination/empty wells":
                            return -12
                        elif row.annotations[1]["value"]=="Artefacts":
                            return -14
                        elif row.annotations[1]["value"]=="Insufficient growth":
                            return -15
                        elif row.annotations[1]["value"]=="Other":
                            return -16
                        else:
                            raise ValueError("unrecognised answer for cannot classify: ",row.annotations[1]["value"])
                    elif row.annotations[1]["value"].isnumeric():
                        return int(row.annotations[1]["value"])
                    else:
                        return(-106)
