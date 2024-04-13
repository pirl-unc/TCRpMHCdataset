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

# Get the directory of the current file (constants.py)
_current_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the absolute path to the 'refs' directory
_refs_dir = os.path.join(_current_dir, '..', 'refs')

HLA_SEQUENCE_MAP = pd.read_csv(os.path.join(_refs_dir, '2field_hla_consensus_seqs.csv'), index_col=0).to_dict()["Full Sequence"]
HLA_PSEUDO_MAP = pd.read_csv(os.path.join(_refs_dir, 'hla_pseudo_seqs.csv'), index_col=0).to_dict()["pseudo-sequence"]
AA_VOCABULARY = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
