#! /usr/bin/env python

import os

from tqdm import tqdm

import pyniverse

class BashTheBugClassifications(pyniverse.Classifications):

    def create_measurements_table(self):

        # self.classifications.drop(['metadata','annotations','subject_data','filename'], axis=1, inplace=True)

        # create a table of measurements, additional measurements (e.g. Vizion or AMyGDA) can be merged in later
        self.measurements=self.classifications[['plate_image','drug','bashthebug_dilution']].groupby(['plate_image','drug']).agg({'bashthebug_dilution':['median','mean','std','min','max','count']})

        # rename the top level of the columns
        self.measurements.columns.set_levels(['bashthebug_dilution'],level=0,inplace=True)


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

        tqdm.pandas(desc='extracting filename ')
        self.classifications['filename']=self.classifications.progress_apply(self._extract_filename2,axis=1)

        # print(self.classifications['filename'])#.str.split('-zooniverse-', 1).str)

        self.classifications['plate_image'], self.classifications['drug']=self.classifications['filename'].str.split('-zooniverse-', 1).str

        tqdm.pandas(desc='calculating dilution')
        self.classifications['bashthebug_dilution']=self.classifications.progress_apply(self._parse_annotation,axis=1).astype(int)

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

    def filter_study(self,study):

        self.classifications=self.classifications.loc[self.classifications["study_id"]==study]

    def filter_readingday(self,reading_day):

        self.classifications=self.classifications.loc[self.classifications["reading_day"]==reading_day]


    def _extract_filename2(self,row):
        try:
            for i in row.subject_data[str(row.subject_ids)]:
                if (".png" in i) or i in ["Filename","Image"]:
                    return(row.subject_data[str(row.subject_ids)][i][:-4])
        except:
            print("Problem parsing "+row.classification_id)
        # print(row)

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
