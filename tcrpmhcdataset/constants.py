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
Keep biological constants here.
"""
import os
import pandas as pd
from importlib import resources
import importlib.resources as pkg_resources

def get_resource_path(filename):
    with pkg_resources.path('tcrpmhcdataset.refs', filename) as path:
        return str(path)

HLA_SEQUENCE_MAP = pd.read_csv(get_resource_path('2field_hla_consensus_seqs.csv'), 
                              index_col=0).to_dict()["Full Sequence"]
HLA_PSEUDO_MAP = pd.read_csv(get_resource_path('hla_pseudo_seqs.csv'), 
                            index_col=0).to_dict()["pseudo-sequence"]
AA_VOCABULARY = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
