#!pip install PyTDC

from tdc.multi_pred import DTI
import pandas as pd

# get current directory
import os

print(os.getcwd())

data = DTI(name="DAVIS")
df = data.get_data()
df.to_csv("local data/Davis/davis_all-KD_SequenceBased_KDregression.csv", index=False)
