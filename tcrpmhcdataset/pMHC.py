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
The purpose of this python3 script is to implement the pMHC dataclass.
"""

from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional, Set
import re
import mhcgnomes
from .constants import *
import tidytcells as tdtc


@dataclass(frozen=True)
class pMHC:
    """
    Define a Meaningful pMHC Class. 

    Args:
       * peptide (str): The peptide sequence
       * hla_allele (str): The HLA allele
       * cognate_tcr (TCR): The cognate TCR
       * reference (str): The reference for this pMHC

    Attributes:
       * peptide (str): The peptide sequence
       * allele (str): The HLA allele
       * tcrs (set): The set of cognate TCRs
       * mhc (str): The MHC sequence
       * pseudo (str): The pseudo MHC sequence
       * references (set): The set of references for this pMHC

    Implements:
       * __post_init__: Initialize the pMHC object
       * __repr__: Return a string representation of the pMHC object for tokenization purposes
       * __str__: Return a string representation of the pMHC object for user interaction
       * __eq__: Check if two pMHC objects are equal
       * __hash__: Return the hash of the pMHC object
       * add_tcr: Add a cognate TCR to the set of cognate TCRs for this pMHC
       * get_tcrs: Get the set of cognate TCRs for this pMHC
       * add_reference: Add a reference to the set of references for this pMHC
       * get_references: Get the set of references for this pMHC
       * hla_allele_parser: Custom HLA allele to standardize the allele and impute canonical HLA from Haplotypes
       * check_mutations: Check the locus of the mutations before making mutations to base sequence
       * apply_mutations: Apply mutations to the base sequence
       * hla_allele2seq: Take a MHCgnomes standardized allele name and return the IMGT, HLAdb sequence
       * hla_allele2pseudo: Take an imperfect allele name and return the NetMHC Pseudo-sequence
    """
    
    # Immutable fields
    peptide: Optional[str] = None
    hla_allele: Optional[str] = None
    cognate_tcr: Optional[object] = None
    reference: Optional[str] = None
    use_pseudo: bool = True
    use_mhc: bool = False
    eager_impute: bool = False
    
    # Mutable fields 
    tcrs: Set[object] = field(default_factory=set, compare=False, hash=False)
    references: Set[str] = field(default_factory=set, compare=False, hash=False)
    
    # Properties to Compute
    @cached_property
    def allele(self):
        # Parse the allele using custom parser
        if self.eager_impute:
            try:
                parsed_allele = mhcgnomes.parse(self.hla_allele)
                assert isinstance(parsed_allele, mhcgnomes.Allele)
            except:
                parsed_allele = mhcgnomes.parse(self.hla_allele_parser(self.hla_allele))
        else:
            parsed_allele = mhcgnomes.parse(self.hla_allele)
        
        if not isinstance(parsed_allele, mhcgnomes.Allele):
            # Try another time to see if it worked
            parsed_allele = mhcgnomes.parse(self.hla_allele_parser(parsed_allele.to_string()))
            if not isinstance(parsed_allele, mhcgnomes.Allele):
                raise ValueError("Could not parse allele")
        
        # Make it a 2-field allele
        parsed_allele.allele_fields = parsed_allele.allele_fields[:2]
        return parsed_allele.to_string()
        
    @cached_property
    def mhc(self):
        # Return the corresponding MHC sequence of the
        # HLA-Allele  w/ mutations included
        return self.hla_allele2seq()
    
    @cached_property
    def pseudo(self):
        # Return the corresponding MHC psuedo-seq of the 
        # HLA-Allele w/ mutations included
        return self.hla_allele2pseudo()
    
    def __post_init__(self):
        """
        Initialize the pMHC object
        """
        input_peptide = self.peptide
        object.__setattr__(self, 'peptide', tdtc.aa.standardize(self.peptide))

        # Initialize mutable sets with cognate TCR and reference, if provided
        if self.cognate_tcr:
            self.tcrs.add(self.cognate_tcr)
        if self.reference:
            self.add_reference(self.reference)

        if self.allele is None or self.peptide is None:
            raise ValueError("Could not complete pMHC creation for peptide: {input_peptide} and allele {self.hla_allele}")
    
    def __repr__(self):
        # Return a string representation of the pMHC object for tokenization purposes
        return f'pMHC(peptide="{self.peptide}", hla_allele="{self.allele}", reference={self.get_references()}, use_pseudo={self.use_pseudo}, use_mhc={self.use_mhc})'
    
    def __str__(self):
        # Return a string representation of the pMHC object for user interaction
        if self.use_mhc:
            return f'{self.peptide}[SEP]{self.mhc}'
        elif self.use_pseudo:
            return f'{self.peptide}[SEP]{self.pseudo}'
        else:
            return self.peptide
    
    def __eq__(self, other):
        return self.peptide == other.peptide and self.mhc == other.mhc
    
    def __hash__(self):
        return hash((self.peptide, self.mhc))

    def add_tcr(self, cognate_tcr):
        """
        Add a cognate TCR to the set of cognate TCRs for this pMHC

        Args:
           * cognate_tcr (TCR): The TCR to add

        Returns:
           * None
        """
        self.tcrs.add(cognate_tcr)
        
    def get_tcrs(self):
        """
        Get the set of cognate TCRs for this pMHC

        Returns:
           * tcrs: The set of cognate TCRs for this pMHC 
        """
        return self.tcrs

    def add_reference(self, reference):
        """
        If another reference is supports this PMHC, add it to the existing set of references.

        Args:
           * reference (str): The reference to add. Can be in any format recognized by the user.

        Returns:
           * None
        """
        if isinstance(reference, str):
            self.references.add(reference)
        elif isinstance(reference, set) or isinstance(reference, list):
            self.references.update(reference)
        else:
            pass

    def get_references(self):
        """
        Get the set of references for this pMHC

        Returns:
           * references: The set of references for this pMHC
        """
        return self.references
        
    def hla_allele_parser(self, hla_string):
        """
        # NOTE: There are weird edge cases that cause failure to impute. 
        Enthusiastic parser. Match the HLA gene and if the serotype information
        is given, cycle through the allele possibilities (1,11) to identify
        the canonical serotype as the first matched one from the available
        sequenced MHCs.

        Args:
           * hla_string (str): The HLA information 

        Returns:
           * parsed_hla_string (str): The standardized HLA allele string
        """
        #Strip whitespace from the input string
        hla_string = hla_string.strip()
        # Add "HLA-" at the beginning of the input string if it is missing
        if not hla_string.startswith("HLA-"):
            hla_string = "HLA-" + hla_string
        # Add an asterisk after the HLA class if it is missing
        if not re.search(r'HLA-[ABC]\*', hla_string):
            hla_string = re.sub(r'HLA-([ABC])', r'HLA-\1*', hla_string)        

        # Regular expression pattern to match HLA allele strings
        pattern = r'^HLA-([ABC])\*(\d+)(?::(\d+))?(?::\d+)?$'
        # Search for matches in the input string
        match = re.search(pattern, hla_string)
        if match:
            # Extract the matched HLA class, allele group, and protein number (if present)
            hla_class = match.group(1)
            allele_group = int(match.group(2))
            for i in range(1,11):
                protein_num = int(match.group(3)) if match.group(3) else i
                if f"HLA-{hla_class}*{allele_group:02d}:{protein_num:02d}" in HLA_SEQUENCE_MAP.keys():
                    return f"HLA-{hla_class}*{allele_group:02d}:{protein_num:02d}"
                 # Return the standardized HLA allele representation
            return f"HLA-{hla_class}*{allele_group:02d}:{protein_num:02d}"
        else:
            # If the input string does not match the pattern, return the original string
            return hla_string

    @staticmethod
    def check_mutations(mutations, og_seq):
        """
        Check the locus of the mutations before making mutations to base sequence
        
        Args: 
           * mutations - Tuple of (position, aa_original, aa_mutant)
           * og_seq - WT consensus sequence [May or may not have linker sequence]
        
        Returns: 
           * True if matches in all loci, false otherwise
        """
        for mutation in mutations:
            pos, orig, mut = mutation
            if og_seq[pos] == orig:
                continue
            else:
                return False
        return True

    @staticmethod
    def apply_mutations(mutations, og_seq):
        """
        Apply mutations to the base sequence.

        Args:
           * mutations - Tuple of (position, aa_original, aa_mutant)
           * og_seq - WT consensus sequence [May or may not have linker sequence]   
        
        Returns:
           * seq - The mutated sequence
        """
        seq = og_seq
        for mutation in mutations:
            pos, orig, mut = mutation
            seq = seq[:pos] + mut + seq[pos + 1 :]
        return seq

    def hla_allele2seq(self):
        """
        Take a MHCgnomes standardized allele name and return the IMGT, HLAdb sequence.
        Capable of Handling mutations through MHCgnomes parser.
        DISCLAIMER: Potential error with C*03:03

        Returns:
           * seq (str): The sequence of the MHC allele
        """
        # To handle mutations or base case, take the root of the allele
        base_allele = str.split(self.allele, " ")[0]
        try:
            seq = HLA_SEQUENCE_MAP[base_allele]
        except KeyError as e:
            print(
                "Could not find MHC reference for {}".format(
                    base_allele
                )
            )
            print("Returning blank sequence.")
            # Returns blank string to be padded out
            return ""

        mhc_allele = mhcgnomes.parse(self.allele)
        # Make any changes to the sequence per mutation
        if len(mhc_allele.mutations) > 0:
            try:
                # Account for zero indexing. Check mutations are aligned before making them
                mutations = [
                    (mut.pos - 1, mut.aa_original, mut.aa_mutant)
                    for mut in mhc_allele.mutations
                ]  # -1 for 0 index
                assert self.check_mutations(mutations, seq)
                seq = self.apply_mutations(mutations, seq)
            except AssertionError:
                try:
                    # Account for signal peptide and 0-index. Check mutation alignment before changing sequence.
                    mutations = [
                        (mut.pos + 24 - 1, mut.aa_original, mut.aa_mutant)
                        for mut in mhc_allele.mutations #for mut in self.allele.mutations
                    ]
                    assert self.check_mutations(mutations, seq)
                    seq = self.apply_mutations(mutations, seq)
                except AssertionError:
                    print(
                        "Could not align mutations to reference for {}".format(
                            self.allele
                        )
                    )
                    print("Defaulting to available reference sequence.")
                    seq = seq
        return seq

    def hla_allele2pseudo(self):
        """
        Take an imperfect allele name and return the NetMHC Pseudo-sequence
        # NOTE: Does not contain values for all the Alleles or handle mutations given predefined pseudo-sequences.

        Returns:
           * pseudo_seq (str): The pseudo-sequence of the MHC allele
        """
        # Assume that mutations do not affect pseudo-sequence
        base_allele = str.split(self.allele, " ")[0]
        try:
            pseudo_seq = HLA_PSEUDO_MAP[base_allele]
        except KeyError:
            # Pseudo-seq reference does not have the 2-field HLA Allele
            print(
                "Could not find PSEUDO reference for {}".format(
                    self.allele
                )
            )
            print("Returning blank sequence.")
            # Returns blank string to be padded out
            pseudo_seq = ""
        return pseudo_seq
    