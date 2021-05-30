import csv
import pandas as pd

score_clf = []
with open("_chemprop_clf.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        score_clf += [float(r[1])]

score_reg = []
with open("_chemprop_reg.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        score_reg += [float(r[1])]

df_ = pd.DataFrame({
    "ActivityClfGraph": score_clf,
    "ActivityRegGraph": score_reg
    }
)

df_.to_csv("data_5.csv", index=False)
