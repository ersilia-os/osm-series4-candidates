# conda activate chemprop

python _5_1_chemprop.py

chemprop_predict --test_path _chemprop.csv \
                 --checkpoint_dir models/chemprop/clf_checkpoints \
                 --preds_path _chemprop_clf.csv

chemprop_predict --test_path _chemprop.csv \
                 --checkpoint_dir models/chemprop/reg_checkpoints \
                 --preds_path _chemprop_reg.csv

python _5_2_chemprop.py
