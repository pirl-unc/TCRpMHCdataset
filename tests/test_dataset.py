###############################################################################
######## Tests for Classes and Functions for TCRpMHCDataset Object ############
###############################################################################

# Import custom classes
from tcrpmhcdataset.dataset import *
import pytest
import os
import pandas as pd
import torch


class MockTokenizer:
    """
    Mock tokenizer for testing collate_fn and get_dataloader.
    Mimics a HuggingFace-style tokenizer interface.
    """
    def __call__(self, text_list, padding=True, truncation=True, max_length=25, return_tensors='pt'):
        """
        Mock tokenization - converts strings to simple integer encodings.
        Returns a dict with 'input_ids' tensor.
        """
        encoded = []
        for text in text_list:
            text = text or ""
            # Simple encoding: convert each character to its ASCII value
            ids = [ord(c) for c in text[:max_length]]
            # Pad to max_length
            ids = ids + [0] * (max_length - len(ids))
            encoded.append(ids)
        return {'input_ids': torch.tensor(encoded)}

curr_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(curr_dir)

class TestTCRpMHCdataset:
    @classmethod
    def setup_class(cls):
        """
        Called before class initialization.
        """
        pass

    @classmethod
    def teardown_class(cls):
        """
        Called after every class initialization.
        """
        pass

    def setup_method(self):
        """
        Called before every method.
        """
        self.tcr2pmhc_dataset = TCRpMHCdataset(
                                source='tcr',
                                target='pmhc',
                                use_mhc=False,
                                use_pseudo=True,
                                use_cdr3=True,
        )
        self.pmhc2tcr_dataset = TCRpMHCdataset(
                                source='pmhc',
                                target='tcr',
                                use_mhc=False,
                                use_pseudo=True,
                                use_cdr3=True,
        )
        self.sample_data_path = os.path.join(curr_dir, 'test_data/sampled_paired_data_cleaned.csv')
        self.sample_data_df = pd.read_csv(self.sample_data_path)

    def teardown_method(self):
        """
        Called after every method.
        """
        del self.tcr2pmhc_dataset
        del self.pmhc2tcr_dataset
        del self.sample_data_df
        del self.sample_data_path

    def test_tcr2pmhc_init(self):
        """
        Test initialization of TCRpMHCdataset object.
        """
        assert self.tcr2pmhc_dataset.source == 'tcr'
        assert self.tcr2pmhc_dataset.target == 'pmhc'
        assert self.tcr2pmhc_dataset.use_mhc == False
        assert self.tcr2pmhc_dataset.use_pseudo == True
        assert self.tcr2pmhc_dataset.use_cdr3 == True
        assert len(self.tcr2pmhc_dataset) == 0
        assert str(self.tcr2pmhc_dataset) == 'TCR:pMHC Dataset of N=0. Mode:tcr -> pmhc.'

    def test_pmhc2tcr_init(self):
        """
        Test initialization of TCRpMHCdataset object.
        """
        assert self.pmhc2tcr_dataset.source == 'pmhc'
        assert self.pmhc2tcr_dataset.target == 'tcr'
        assert self.pmhc2tcr_dataset.use_mhc == False
        assert self.pmhc2tcr_dataset.use_pseudo == True
        assert self.pmhc2tcr_dataset.use_cdr3 == True
        assert len(self.pmhc2tcr_dataset) == 0
        assert str(self.pmhc2tcr_dataset) == 'TCR:pMHC Dataset of N=0. Mode:pmhc -> tcr.'

    def test_get_srclist(self):
        """
        Test getting the source list of the dataset.
        """
        self.tcr2pmhc_dataset.load_data_from_df(self.sample_data_df)
        self.pmhc2tcr_dataset.load_data_from_df(self.sample_data_df)
        for src in self.tcr2pmhc_dataset.get_srclist():
            assert isinstance(src, TCR)
        for src in self.pmhc2tcr_dataset.get_srclist():
            assert isinstance(src, pMHC)

    def test_get_trglist(self):
        # Test getting the source list of the dataset.
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        for trg in self.tcr2pmhc_dataset.get_trglist():
            assert isinstance(trg, pMHC)
        for trg in self.pmhc2tcr_dataset.get_trglist():
            assert isinstance(trg, TCR)

    def test_to_dict(self):
        # Test the sequence to sequence mapping dict generated from the dataset
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        
        tcr2pmhc_dict = self.tcr2pmhc_dataset.to_dict()
        pmhc2tcr_dict = self.pmhc2tcr_dataset.to_dict()
        

        # Test TCR to pMHC mapping
        for src, trg_set in tcr2pmhc_dict.items():
            # Check TCR expectations
            assert isinstance(src, TCR)
            assert '[SEP]' not in str(src)
            assert str(src).startswith('C')

            # Check pMHC expectations
            assert isinstance(trg_set, list)
            for trg in trg_set:
                assert isinstance(trg, str)
                assert '[SEP]' in trg

    
        # Test pMHC to TCR mapping
        for src, trg_set in pmhc2tcr_dict.items():
            # Check pMHC expectations
            assert isinstance(src, pMHC)
            assert '[SEP]' in str(src)

            # Check pMHC expectations
            assert isinstance(trg_set, list)
            for trg in trg_set:
                assert isinstance(trg, str)
                assert '[SEP]' not in trg
                assert trg.startswith('C')
    
    def test_load_data_from_file(self):
        """
        Test loading of data from a File by first making it into a suitable DF.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        assert len(self.tcr2pmhc_dataset) == 6833
        assert str(self.tcr2pmhc_dataset) == 'TCR:pMHC Dataset of N=6833. Mode:tcr -> pmhc.'
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        assert len(self.pmhc2tcr_dataset) == 6833
        assert str(self.pmhc2tcr_dataset) == 'TCR:pMHC Dataset of N=6833. Mode:pmhc -> tcr.'

    def test_load_data_from_df(self):
        """
        Test loading of data from a suitable DF.
        """
        self.tcr2pmhc_dataset.load_data_from_df(self.sample_data_df)
        assert len(self.tcr2pmhc_dataset) == 6833
        assert str(self.tcr2pmhc_dataset) == 'TCR:pMHC Dataset of N=6833. Mode:tcr -> pmhc.'
        self.pmhc2tcr_dataset.load_data_from_df(self.sample_data_df)
        assert len(self.pmhc2tcr_dataset) == 6833
        assert str(self.pmhc2tcr_dataset) == 'TCR:pMHC Dataset of N=6833. Mode:pmhc -> tcr.'

    def test_to_df(self):
        """
        Test conversion of TCRpMHCdataset to pandas dataframe.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
    
        assert isinstance(self.tcr2pmhc_dataset.to_df(), pd.DataFrame)
        assert len(self.tcr2pmhc_dataset.to_df()) == len(self.tcr2pmhc_dataset)
        assert isinstance(self.tcr2pmhc_dataset, TCRpMHCdataset)

        tcr_cdr3bs = [tcr.cdr3b for tcr in self.tcr2pmhc_dataset.tcrs]
        tcr_trbvs = [tcr.trbv for tcr in self.tcr2pmhc_dataset.tcrs]
        tcr_trbjs = [tcr.trbj for tcr in self.tcr2pmhc_dataset.tcrs]

        pmhc_peptides = [pmhc.peptide for pmhc in self.tcr2pmhc_dataset.pMHCs]
        pmhc_alleles = [pmhc.allele for pmhc in self.tcr2pmhc_dataset.pMHCs]

        # Assert that all TCR and pMHC sequences are in the dataset and in the same spot
        assert list(self.tcr2pmhc_dataset.to_df()['CDR3b']) == tcr_cdr3bs
        assert list(self.tcr2pmhc_dataset.to_df()['Epitope']) == pmhc_peptides
        assert list(self.tcr2pmhc_dataset.to_df()['Allele']) == pmhc_alleles
        assert list(self.tcr2pmhc_dataset.to_df()['TRBV']) == tcr_trbvs
        assert list(self.tcr2pmhc_dataset.to_df()['TRBJ']) == tcr_trbjs

    def test_to_csv(self):
        """
        Test conversion of TCRpMHCdataset to csv file.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        self.tcr2pmhc_dataset.to_csv('test_data.csv')
        assert os.path.exists('test_data.csv')
        os.remove('test_data.csv')

    def test_getitem(self): 
        """
        Test getting an item from the dataset.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        assert isinstance(self.tcr2pmhc_dataset[0], tuple)
        assert isinstance(self.pmhc2tcr_dataset[0], tuple)
        assert isinstance(self.tcr2pmhc_dataset[0][0], TCR)
        assert isinstance(self.tcr2pmhc_dataset[0][1], pMHC)
        assert isinstance(self.pmhc2tcr_dataset[0][0], pMHC)
        assert isinstance(self.pmhc2tcr_dataset[0][1], TCR)

    def test_unbalanced_split(self):
        """
        Test splitting of dataset into train and test sets.
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        train_dset, test_dset = self.pmhc2tcr_dataset.split(test_size=0.2, balance_on_allele=False)

        # Check that the split gives the correct number of samples each
        assert len(train_dset)/(len(train_dset) + len(test_dset)) == pytest.approx(0.8, .05)
        assert len(test_dset)/(len(train_dset) + len(test_dset)) == pytest.approx(0.2, .05)

    def test_balanced_split(self):
        """
        Test splitting of dataset into train and test sets.
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        train_dset, test_dset = self.pmhc2tcr_dataset.split(test_size=0.2, balance_on_allele=True)

        train_df = train_dset.to_df()
        test_df = test_dset.to_df()

        # Check that the split gives the correct number of samples each
        assert len(train_df)/(len(train_df) + len(test_df)) == pytest.approx(0.8, .1)
        assert len(test_df)/(len(train_df) + len(test_df)) == pytest.approx(0.2, .1)

        # Check that the splits are balanced
        train_alleles_df = train_df.groupby('Allele').count()
        train_alleles_freq = train_alleles_df['CDR3b']/train_alleles_df['CDR3b'].sum()
        test_alleles_df = test_df.groupby('Allele').count()
        test_alleles_freq = test_alleles_df['CDR3b']/test_alleles_df['CDR3b'].sum()

        # Since all singletons go to the train_df 
        # check that all of the test_df Alleles are in the train_df
        assert set(test_alleles_df.index).issubset(set(train_alleles_df.index))
        for allele in test_alleles_df.index:
            assert test_alleles_freq[allele] == pytest.approx(train_alleles_freq[allele], abs=.1)
        

    def test_balanced_split_on_epitope(self):
        """
        Test splitting the dataset into train and test sets based on epitope (should not be repeated in train/test).
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        train_dset, test_dset = self.pmhc2tcr_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['Epitope'])

        train_df = train_dset.to_df()
        test_df = test_dset.to_df()

        # Check that the split gives the correct number of samples each
        assert len(train_df)/(len(train_df) + len(test_df)) == pytest.approx(0.8, abs=.1)
        assert len(test_df)/(len(train_df) + len(test_df)) == pytest.approx(0.2, abs=.1)

        # Check that the splits are balanced
        train_alleles_df = train_df.groupby('Allele').count()
        train_alleles_freq = train_alleles_df['CDR3b']/train_alleles_df['CDR3b'].sum()
        test_alleles_df = test_df.groupby('Allele').count()
        test_alleles_freq = test_alleles_df['CDR3b']/test_alleles_df['CDR3b'].sum()

        # Since all singletons go to the train_df 
        # check that all of the test_df Alleles are in the train_df
        for allele in test_alleles_df.index:
            if allele in train_alleles_df.index:
                assert test_alleles_freq[allele] == pytest.approx(train_alleles_freq[allele], abs=.1)

        assert len(set(train_df['Epitope']).intersection(set(test_df['Epitope']))) == 0

    def test_balanced_split_on_pmhc(self):
        """
        Test splitting the dataset into train and test sets based on pmhc (should not be repeated in train/test).
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        train_dset, test_dset = self.pmhc2tcr_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['Epitope', 'Allele'])

        train_df = train_dset.to_df()
        test_df = test_dset.to_df()

        # Check that the split gives the correct number of samples each
        assert len(train_df)/(len(train_df) + len(test_df)) == pytest.approx(0.8, abs=.1)
        assert len(test_df)/(len(train_df) + len(test_df)) == pytest.approx(0.2, abs=.1)

        # Check that the splits are balanced
        train_alleles_df = train_df.groupby('Allele').count()
        train_alleles_freq = train_alleles_df['CDR3b']/train_alleles_df['CDR3b'].sum()
        test_alleles_df = test_df.groupby('Allele').count()
        test_alleles_freq = test_alleles_df['CDR3b']/test_alleles_df['CDR3b'].sum()

        # Since all singletons go to the train_df 
        # check that all of the test_df Alleles are in the train_df
        assert set(test_alleles_df.index).issubset(set(train_alleles_df.index))
        for allele in test_alleles_df.index:
            if allele in train_alleles_df.index:
                assert test_alleles_freq[allele] == pytest.approx(train_alleles_freq[allele], abs=.1)

        train_df['group'] = train_df['Epitope'] + '_' + train_df['Allele']
        test_df['group'] = test_df['Epitope'] + '_' + test_df['Allele']
        assert len(set(train_df['group']).intersection(set(test_df['group']))) == 0

    def test_balanced_split_on_tcr(self):
        """
        Test splitting the dataset into train and test sets based on tcr (should not be repeated in train/test).
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        train_dset, test_dset = self.pmhc2tcr_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['CDR3b', 'TRBV', 'TRBJ'])

        train_df = train_dset.to_df()
        test_df = test_dset.to_df()

        # Check that the split gives the correct number of samples each
        assert len(train_df)/(len(train_df) + len(test_df)) == pytest.approx(0.8, abs=.1)
        assert len(test_df)/(len(train_df) + len(test_df)) == pytest.approx(0.2, abs=.1)

        # Check that the splits are balanced
        train_alleles_df = train_df.groupby('Allele').count()
        train_alleles_freq = train_alleles_df['CDR3b']/train_alleles_df['CDR3b'].sum()
        test_alleles_df = test_df.groupby('Allele').count()
        test_alleles_freq = test_alleles_df['CDR3b']/test_alleles_df['CDR3b'].sum()

        # Since all singletons go to the train_df 
        # check that all of the test_df Alleles are in the train_df
        for allele in test_alleles_df.index:
            if allele in train_alleles_df.index:
                assert test_alleles_freq[allele] == pytest.approx(train_alleles_freq[allele], abs=.1)

        train_df['group'] = train_df['CDR3b'] + '_' + train_df['TRBV'] + '_' + train_df['TRBJ']
        test_df['group'] = test_df['CDR3b'] + '_' + test_df['TRBV'] + '_' + test_df['TRBJ']

        assert len(set(train_df['group']).intersection(set(test_df['group']))) == 0

    def test_collate_fn_tcr2pmhc(self):
        """
        Test the collate function for tcr -> pmhc mode.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        tokenizer = MockTokenizer()
        collate_fn = self.tcr2pmhc_dataset.collate_fn(tokenizer)

        # Get a small batch of items
        batch = [self.tcr2pmhc_dataset[i] for i in range(4)]
        result = collate_fn(batch)

        # Check that result has expected keys
        assert 'input_ids' in result
        assert 'labels' in result

        # Check tensor shapes
        assert result['input_ids'].shape[0] == 4  # batch size
        assert result['labels'].shape[0] == 4  # batch size

        # In tcr2pmhc mode, input_ids should be TCRs and labels should be pMHCs
        assert result['input_ids'].shape[1] == 25  # TCR max_length
        assert result['labels'].shape[1] == 52  # pMHC max_length

    def test_collate_fn_pmhc2tcr(self):
        """
        Test the collate function for pmhc -> tcr mode.
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        tokenizer = MockTokenizer()
        collate_fn = self.pmhc2tcr_dataset.collate_fn(tokenizer)

        # Get a small batch of items
        batch = [self.pmhc2tcr_dataset[i] for i in range(4)]
        result = collate_fn(batch)

        # Check that result has expected keys
        assert 'input_ids' in result
        assert 'labels' in result

        # Check tensor shapes
        assert result['input_ids'].shape[0] == 4  # batch size
        assert result['labels'].shape[0] == 4  # batch size

        # In pmhc2tcr mode, input_ids should be pMHCs and labels should be TCRs
        assert result['input_ids'].shape[1] == 52  # pMHC max_length
        assert result['labels'].shape[1] == 25  # TCR max_length

    def test_get_dataloader_tcr2pmhc(self):
        """
        Test the get_dataloader method for tcr -> pmhc mode.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        tokenizer = MockTokenizer()
        dataloader = self.tcr2pmhc_dataset.get_dataloader(tokenizer, batch_size=8, shuffle=False)

        # Check that dataloader is correct type
        assert isinstance(dataloader, torch.utils.data.DataLoader)

        # Get first batch and verify structure
        batch = next(iter(dataloader))
        assert 'input_ids' in batch
        assert 'labels' in batch
        assert batch['input_ids'].shape[0] == 8
        assert batch['labels'].shape[0] == 8

    def test_get_dataloader_pmhc2tcr(self):
        """
        Test the get_dataloader method for pmhc -> tcr mode.
        """
        self.pmhc2tcr_dataset.load_data_from_file(self.sample_data_path)
        tokenizer = MockTokenizer()
        dataloader = self.pmhc2tcr_dataset.get_dataloader(tokenizer, batch_size=8, shuffle=False)

        # Check that dataloader is correct type
        assert isinstance(dataloader, torch.utils.data.DataLoader)

        # Get first batch and verify structure
        batch = next(iter(dataloader))
        assert 'input_ids' in batch
        assert 'labels' in batch
        assert batch['input_ids'].shape[0] == 8
        assert batch['labels'].shape[0] == 8

    def test_dataloader_iteration(self):
        """
        Test that the dataloader can iterate through multiple batches.
        """
        self.tcr2pmhc_dataset.load_data_from_file(self.sample_data_path)
        tokenizer = MockTokenizer()
        batch_size = 32
        dataloader = self.tcr2pmhc_dataset.get_dataloader(tokenizer, batch_size=batch_size, shuffle=False)

        # Count batches and total samples
        total_samples = 0
        batch_count = 0
        for batch in dataloader:
            batch_count += 1
            total_samples += batch['input_ids'].shape[0]

        # Verify we processed all samples
        assert total_samples == len(self.tcr2pmhc_dataset)
        expected_batches = (len(self.tcr2pmhc_dataset) + batch_size - 1) // batch_size
        assert batch_count == expected_batches
