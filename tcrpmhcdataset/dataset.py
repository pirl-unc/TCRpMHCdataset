# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement the TCRpMHCdataset class.
"""

import pandas as pd
import numpy as np
from .constants import *
from .TCR import TCR
from .pMHC import pMHC
import warnings
from sklearn.model_selection import train_test_split

class TCRpMHCdataset:
    """
    Main class for the TCRpMHCDataset package. This class is designed to take tabular paired data as input and return 
    a cohesive dataset that is designed to capture the many to many nature of TCR and pMHC cross-reactivity. Accepts 
    either TCR -> multiple pMHC mapping or pMHC -> multiple TCR mapping. Parses data from a table into lists of TCR and pMHC 
    objects which can be indexed and called during training/eval. The dataset can also be split into stratified train/test sets. 

    
    Args: 
       * source (str): The source of the dataset. Either 'tcr' or 'pmhc'. 
       * target (str): The target of the dataset. Either 'tcr' or 'pmhc'. 
       * use_mhc (bool): Whether to use the MHC sequence or the pMHC sequence. 
       * use_pseudo (bool): Whether to use the pseudo MHC sequence or the full MHC sequence. 
       * use_cdr3 (bool): Whether to use the CDR3 sequence or the full TCR sequence. 
    
        
    Attributes: 
       * source (str): The source of the dataset. Either 'tcr' or 'pmhc'.
       * target (str): The target of the dataset. Either 'tcr' or 'pmhc'.
       * tcrs (list): The list of TCRs in the dataset.
       * pMHCs (list): The list of pMHCs in the dataset.
       * use_mhc (bool): Whether to use the MHC sequence or the pMHC sequence.
       * use_pseudo (bool): Whether to use the pseudo MHC sequence or the full MHC sequence.
       * use_cdr3 (bool): Whether to use the CDR3 sequence or the full TCR sequence.
       * use_both_chains (bool): Whether to use both alpha and beta chains of the TCR.

    Implements:
       * __len__ (self): Return the number of TCR:pMHC pairs in the dataset. Can be accessed using len(dataset).
       * __repr__ (self): Return a string representation of the dataset. Can be accessed using repr(dataset).
       * __str__ (self): Return a more user friendly string representation of the dataset. Can be accessed using str(dataset).
       * __getitem__ (self, idx): Return a tuple of (source object, target object) for the given index. Is either (TCR, PMHC)
        or (PMHC, TCR) depending on the source and target attributes. Can be accessed using dataset[idx].
       * get_srclist (self): Return the list of source objects.
       * get_trglist (self): Return the list of target objects.
       * load_data_from_file (self, path_to_csv): Load the data from a csv file
       * load_data_from_df (self, df): Load the data from a dataframe
       * split (self, test_size=0.2, balance_on_allele=True, split_on=None, random_seed=42): Split the dataset into train, and test data.
       * to_dict (self, stringify_input=False, stringify_output=True): Return a de-dpuplicated dictionary representation of the parallel dataset.
       * to_df (self): Return a dataframe representation of the dataset instance.
       * to_csv (self, path_to_csv): Write the dataset to a csv file.
    """

    def __init__(self, source, target, use_mhc=False, use_pseudo=True, use_cdr3=True, use_both_chains=False):
        assert source in ['tcr', 'pmhc']
        assert target in ['tcr', 'pmhc']
        self.source = source
        self.target = target
        self.tcrs = []
        self.pMHCs = []
        self._tcr_dict = dict()
        self._pmhc_dict = dict()
        self.use_mhc = use_mhc
        self.use_pseudo = use_pseudo if not use_mhc else False
        self.use_cdr3 = use_cdr3
        self.use_both_chains = use_both_chains
        
    def __len__(self):
        """Return the number of TCR:pMHC pairs in the dataset."""
        assert len(self.pMHCs) == len(self.tcrs)
        return len(self.pMHCs)
    
    def __repr__(self):
        """Return a string representation of the dataset."""
        return f'TCRpMHCdataset(source="{self.source}", target="{self.target}",use_mhc={self.use_mhc},use_pseudo={self.use_pseudo}, use_cdr3={self.use_cdr3})'
    
    def __str__(self):
        """Return a more user friendly string representation of the dataset."""
        return f'TCR:pMHC Dataset of N={self.__len__()}. Mode:{self.source} -> {self.target}.'
    
    def __getitem__(self, idx):
        """Return a tuple of (source object, target object) for the given index. Is either (TCR, PMHC)
        or (PMHC, TCR) depending on the source and target attributes.
        """
        tcr = self.tcrs[idx] 
        pmhc = self.pMHCs[idx]
        if self.source == 'pmhc':
            return pmhc, tcr
        else:    
            return tcr, pmhc
        
    def get_srclist(self):
        """Return the list of source objects."""
        return self.tcrs if self.source == 'tcr' else self.pMHCs
    
    def get_trglist(self):
        """Return the list of target objects."""
        return self.pMHCs if self.source == 'tcr' else self.tcrs
    
    def load_data_from_file(self, path_to_csv, verbose=False):
        """
        Load the data from a csv file with the following required columns:

            1. 'CDR3b'
            2. 'TRBV'
            3. 'TRBJ'
            4. 'Epitope'
            5. 'Allele'
            6. 'Reference'
        
        Args:
           * path_to_csv (str): The path to the csv file with the following columns:

                1. 'CDR3a': The CDR3a sequence in capital single letter Amino Acid Code format (str, optional)
              * 2. 'CDR3b': The CDR3b sequence in capital single letter Amino Acid Code format (str, required)
                3. 'TRAV': The TRAV gene in IMGT format (str, optional)
              * 4. 'TRBV': The TRBV gene in IMGT format (str, required)
                5. 'TRAJ': The TRAJ gene in IMGT format (str, optional)
              * 6. 'TRBJ': The TRBJ gene in IMGT format (str, required)
                7. 'TRAD': The TRAD gene in IMGT format (str, optional)
                8. 'TRBD': The TRBD gene in IMGT format (str, optional)
                9. 'TRA_stitched': The full TRA sequence in capital single letter Amino Acid Code format (str, optional)
                10. 'TRB_stitched': The full TRB sequence in capital single letter Amino Acid Code format (str, optional [can be imputed])
              * 11. 'Epitope': The peptide sequence in capital single letter Amino Acid Code format (str, required)
              * 12. 'Allele': The HLA allele in IMGT format (str, required)
                13. 'Pseudo': The pseudo MHC sequence in capital single letter Amino Acid Code format (str, optional [can be imputed])
                14. 'MHC': The full MHC sequence in capital single letter Amino Acid Code format (str, optional [can be imputed])
              * 15. 'Reference': The reference for the data point (str, required)
        
        Raises: 
           * FileNotFoundError: If the file is not found. 
           * Warnings if specific instances were unable to be loaded.
        
        Returns:
           * None: This function does not return anything.

        """
        try:
            df = pd.read_csv(path_to_csv)
        except FileNotFoundError:
            print(f'File not found: {path_to_csv}')
            return
        
        assert 'CDR3b' in df.columns
        assert 'TRBV' in df.columns
        assert 'TRBJ' in df.columns
        assert 'Epitope' in df.columns
        assert 'Allele' in df.columns
        assert 'Reference' in df.columns
        
        if 'CDR3a' not in df.columns:
            df['CDR3a'] = ''
        if 'TRAV' not in df.columns:
            df['TRAV'] = ''
        if 'TRAJ' not in df.columns:
            df['TRAJ'] = ''
        if 'TRAD' not in df.columns:
            df['TRAD'] = ''
        if 'TRBD' not in df.columns:
            df['TRBD'] = ''
        if 'TRA_stitched' not in df.columns:
            df['TRA_stitched'] = ''
        if 'TRB_stitched' not in df.columns:
            df['TRB_stitched'] = ''
        if 'Pseudo' not in df.columns:
            df['Pseudo'] = ''
        if 'MHC' not in df.columns:
            df['MHC'] = ''
            
        self.load_data_from_df(df, verbose=verbose)

    def load_data_from_df(self, df, verbose=False):
        """
        Load the data from a dataframe with the following required columns:

            1. 'CDR3a'
            2. 'CDR3b'
            3. 'TRAV'
            4. 'TRBV'
            5. 'TRAJ'
            6. 'TRBJ'
            7. 'TRAD'
            8. 'TRBD'
            9. 'TRA_stitched'
            10. 'TRB_stitched'
            11. 'Epitope'
            12. 'Allele'
            13. 'Pseudo'
            14. 'MHC'
            15. 'Reference'

        Returns:
           * None: This function does not return anything.
        """
        og_nrows = len(df)
        df = df.replace({np.nan: None})
        for index, row in df.iterrows():
            try:
                ### 1. Create the TCR and pMHC objects
                use_cdr3a = True if (self.use_cdr3 and self.use_both_chains) else False
                use_trb = True if not self.use_cdr3 else False
                use_tra = True if (~self.use_cdr3 and self.use_both_chains) else False
                tcr_i = TCR(cdr3a=(row['CDR3a'] if isinstance(row['CDR3a'], str) else None), 
                            cdr3b=row['CDR3b'], 
                            trav=row['TRAV'], trbv=row['TRBV'], 
                            traj=row['TRAJ'], trbj=row['TRBJ'],
                            trad=row['TRAD'], trbd=row['TRBD'], 
                            tcra_full=row['TRA_stitched'], tcrb_full=row['TRB_stitched'],
                            reference=row['Reference'], use_cdr3b=self.use_cdr3, use_cdr3a=use_cdr3a, 
                            use_trb=use_trb, use_tra=use_tra)
            except:
                if verbose:
                    warnings.warn(f'Error loading row {index} TCR(cdr3a={row["CDR3a"]}, cdr3b={row["CDR3b"]},'
                              f'trbv={row["TRBV"]},trbj={row["TRBJ"]}). Skipping...', RuntimeWarning)
                continue
            try:
                pMHC_i = pMHC(peptide=row['Epitope'], hla_allele=row['Allele'], reference=row['Reference'], use_pseudo=self.use_pseudo, use_mhc=self.use_mhc, eager_impute=True)
            except:
                if verbose:
                    warnings.warn(f'Error loading row {index} pMHC(peptide={row["Epitope"]}, allele={row["Allele"]}). Skipping...', RuntimeWarning)
                continue
            
            ### 2. Hash the TCR and pMHC objects to get unique keys
            tcr_key = hash(tcr_i)
            pMHC_key = hash(pMHC_i)

            ### 3. If tcr exists then grab the existing tcr object and add the new information to it
            if tcr_key in self._tcr_dict.keys():
                tcr_i = self._tcr_dict[tcr_key]

            # Add reference and cognate pMHC information to that TCR (assumes no duplicates of paired data)
            tcr_i.add_reference(row['Reference'])
            tcr_i.add_pMHC(pMHC_i)
            # Add the updated version back to the dictionary
            self._tcr_dict[tcr_key] = tcr_i

            ### 4. If pmhc exists then grab the existing pmhc object and add the new information to it
            if pMHC_key in self._pmhc_dict.keys():
                pMHC_i = self._pmhc_dict[pMHC_key]
            # Add reference and cognate TCR
            pMHC_i.add_reference(row['Reference'])
            pMHC_i.add_tcr(tcr_i)
            # Add the updated version to the dictionary
            self._pmhc_dict[pMHC_key] = pMHC_i
            
            # Add TCR and PMHC to list **Updates the previous objects in the list thanks to pythons pointers**
            self.tcrs.append(tcr_i)
            self.pMHCs.append(pMHC_i)
            
                
        print(f'Loaded {len(self)} TCR:pMHC pairs from {og_nrows} rows of data.')

    def split(self, test_size=0.2, balance_on_allele=True, split_on=None, random_seed=42):
        """
        Split the dataset into train, and test data. The split is stratified by allele so that the allele distributions of the train and test sets 
        are approximately equal. The train_test_split also contains functionality of ensuring that epitopes and/or TCRs are held out from the training set to 
        assess the generalization capacity of the model.
        
        Args:
           * test_size (float): The proportion of the dataset to include in the test set.
           * balance_on_allele (bool): Whether to balance the train/test split on allele.
           * split_on (list): The column(s) to split on, ensures combinations of instances from these columns occur in both train and test.  
                            *    ['Epitope'] ensures that no epitope is shared between train and test.  
                            *    ['Epitope', 'Allele'] ensures that no epitope::allele combination is shared between train and test.  
                            *    ['CDR3a', 'CDR3b'] ensures that no CDR3a::CDR3b combination is shared between train and test.  
           * random_seed (int): The random seed to use for the train/test split.

        Returns:
           * train_dataset (TCRpMHCdataset): The training dataset.
           * test_dataset (TCRpMHCdataset): The testing dataset.
        """
        assert test_size > 0 and test_size < 1, 'Test size must be between 0 and 1.'
        np.random.seed(random_seed)  # Consistent random seed usage
        data_df = self.to_df()
        train_size = 1 - test_size
        # Handle alleles
        allele_counts = data_df['Allele'].value_counts()
        data_df_singletons = data_df[data_df['Allele'].isin(allele_counts[allele_counts == 1].index)]
        data_df_no_singletons = data_df[~data_df['Allele'].isin(allele_counts[allele_counts == 1].index)]

        # Create split columns if needed
        if split_on:
            data_df_singletons['split_col'] = data_df_singletons[split_on].astype(str).agg('::'.join, axis=1)
            data_df_no_singletons['split_col'] = data_df_no_singletons[split_on].astype(str).agg('::'.join, axis=1)
            unique_split_vals = data_df_no_singletons['split_col'].unique()
            train_vals = np.random.choice(unique_split_vals, size=int(train_size * len(unique_split_vals)), replace=False)
            train_df = data_df_no_singletons[data_df_no_singletons['split_col'].isin(train_vals)]
            test_df = data_df_no_singletons[~data_df_no_singletons['split_col'].isin(train_vals)]
        else:
            train_df, test_df = train_test_split(data_df_no_singletons, test_size=test_size, stratify=data_df_no_singletons['Allele'] if balance_on_allele else None, random_state=random_seed)

        # Adjust for singletons
        if not data_df_singletons.empty:
            train_df = pd.concat([train_df, data_df_singletons], ignore_index=True) # Something is happening here

        # Post Hoc Ensure that Split-On Instances Do Not Co-Occur
        if split_on:
            #test_df = test_df[~test_df['split_col'].isin(train_df['split_col'].unique())]
            #train_df = pd.concat([train_df, test_df[test_df['split_col'].isin(train_df['split_col'].unique())]])
            test_df = pd.concat([test_df, train_df[train_df['split_col'].isin(test_df['split_col'].unique())]])
            train_df = train_df[~train_df['split_col'].isin(test_df['split_col'].unique())]
            
            # Check that the split was successful and split-on instances do not co-occur
            assert len(set(test_df['split_col']).intersection(set(train_df['split_col']))) == 0
            # Get a boolean DataFrame indicating NaN values
            nan_df = train_df.isna()
            # Check if any NaN values exist in each row
            row_has_nan = nan_df.any(axis=1)
            # Filter the original DataFrame to get rows with NaN values
            rows_with_nan = train_df[row_has_nan][['Epitope', 'split_col']]

            # Drop the split column
            train_df = train_df.drop(columns=['split_col']) if 'split_col' in train_df.columns else train_df
            test_df = test_df.drop(columns=['split_col']) if 'split_col' in test_df.columns else test_df

        # Load the data back train and test dataset objects
        train_dataset = TCRpMHCdataset(source=self.source, target=self.target, use_mhc=self.use_mhc, use_pseudo=self.use_pseudo, use_cdr3=self.use_cdr3)
        train_dataset.load_data_from_df(train_df, verbose=False)

        test_dataset = TCRpMHCdataset(source=self.source, target=self.target, use_mhc=self.use_mhc, use_pseudo=self.use_pseudo, use_cdr3=self.use_cdr3)
        test_dataset.load_data_from_df(test_df, verbose=False)

        return train_dataset, test_dataset

    def to_dict(self, stringify_input=False, stringify_output=True):    
        """
        Return a de-dpuplicated dictionary representation of the parallel dataset. 
        Keys are the source objects or their string representations, with string representations
        having more condensing of the data (by merging collisions). 
        
        Args:
           * stringify_input (bool): Whether to convert the input to a string representation using __str__.
           * stringify_output (bool): Whether to convert the output to a string representation using __str__.

        Raises:
           * None

        Returns:
           * data_dict (dict): The dictionary representation of the dataset with 
                                k,v pairs of some combination of source, target and
                                the repr function [{TCR: pMHC} repr(TCR):repr(pMHC)].
        """
        data_dict = dict()
        src_dict = self._tcr_dict if self.source == 'tcr' else self._pmhc_dict
        
        for src in src_dict.values():
            ref_trgs = src.get_pMHCs() if self.source == 'tcr' else src.get_tcrs()
            
            key = src if not stringify_input else str(src)
            values = [str(trg) if stringify_output else trg for trg in ref_trgs]
            
            # Update the dictionary
            if key in data_dict:
                data_dict[key].update(values)
            else:
                data_dict[key] = set(values)

        # Covner the sets to lists
        for k, v in data_dict.items():
            data_dict[k] = list(v)

        return data_dict
    
    def to_df(self):
        """
        Return a dataframe representation of the dataset instance.  
        
        Args:
           * None

        Returns:
           * df (pd.DataFrame): The dataframe representation of the dataset where each row is a unique TCR:pMHC pair and the reference is the  
            concatenation of the list of references for that pair. Multiple references are separated by a semicolon.
        """
        df = pd.DataFrame()
        df['CDR3a'] = [tcr.cdr3a if tcr.cdr3a is not None else '' for tcr in self.tcrs]
        df['CDR3b'] = [tcr.cdr3b if tcr.cdr3b is not None else '' for tcr in self.tcrs]
        df['TRAV'] = [tcr.trav if tcr.trav is not None else '' for tcr in self.tcrs]
        df['TRBV'] = [tcr.trbv if tcr.trbv is not None else '' for tcr in self.tcrs]
        df['TRAJ'] = [tcr.traj if tcr.traj is not None else '' for tcr in self.tcrs]
        df['TRBJ'] = [tcr.trbj if tcr.trbj is not None else '' for tcr in self.tcrs]
        df['TRAD'] = [tcr.trad if tcr.trad is not None else '' for tcr in self.tcrs]
        df['TRBD'] = [tcr.trbd if tcr.trbd is not None else '' for tcr in self.tcrs]
        df['TRA_stitched'] = [tcr.tcra_full if tcr.tcra_full is not None else '' for tcr in self.tcrs]
        df['TRB_stitched'] = [tcr.tcrb_full if tcr.tcrb_full is not None else '' for tcr in self.tcrs]
        df['Epitope'] = [pmhc.peptide if pmhc.peptide is not None else '' for pmhc in self.pMHCs]
        df['Allele'] = [pmhc.allele if pmhc.allele is not None else '' for pmhc in self.pMHCs]
        df['Pseudo'] = [pmhc.pseudo if pmhc.pseudo is not None else '' for pmhc in self.pMHCs]
        df['MHC'] = [pmhc.mhc for pmhc in self.pMHCs]
        df['Reference'] = [pmhc.get_references() for pmhc in self.pMHCs]
        return df
    
    def to_csv(self, path_to_csv):
        """
        Write the dataset to a csv file.
        
        Args:
           * path_to_csv (str): The path to the csv file to write to.


        Returns:
           * None: This function does not return anything.
        """
        df = self.to_df()
        df.to_csv(path_to_csv, index=False)

