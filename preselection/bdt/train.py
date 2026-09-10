import numpy as np

import xgboost as xgb

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt

import ROOT as r

import awkward as ak

data = ak.from_parquet("../data/sig.parquet")

X = np.column_stack([ak.flatten(data[col]).to_numpy() for col in data.fields[:12]])
y = ak.flatten(data["labels"]).to_numpy()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

params = {
    "objective": "binary:logistic",
    "device": "cuda:1",
    "tree_method": "hist",
    "eval_metric": "auc",
    "scale_pos_weight": sum(y_train == 0) / sum(y_train == 1),
    "n_estimators": 500,
    "early_stopping_rounds": 100,
    "learning_rate": 0.03,
    "max_depth": 4,
    "min_child_weight": 6,
    "gamma": 0.5,
    "subsample": 0.85,
    "colsample_bytree": 0.9,
    "colsample_bylevel": 0.9,
    "max_delta_step": 0,
    "lambda": 5.0,
    "alpha": 0.5,
}

bdt = xgb.XGBClassifier(**params)

bdt.fit(X_train, y_train, eval_set=[(X_test, y_test)], verbose=True)

# move bdt to cpu
bdt.set_params(device="cpu")

y_pred_test = bdt.predict_proba(X_test)[:, 1]
y_pred_train = bdt.predict_proba(X_train)[:, 1]

fig, ax = plt.subplots(1,2)
fig.set_size_inches(12, 6)

ax[0].hist(y_pred_test[y_test == 0], bins=20, histtype="step", color="b", label="Non VBS Jets", density=True)
ax[0].hist(y_pred_test[y_test == 1], bins=20, histtype="step", color="r", label="VBS Jets", density=True)
ax[0].set_title("Test Set")
ax[0].legend()

ax[1].hist(y_pred_train[y_train == 0], bins=20, histtype="step", color="b", label="Non VBS Jets", density=True)
ax[1].hist(y_pred_train[y_train == 1], bins=20, histtype="step", color="r", label="VBS Jets", density=True)
ax[1].set_title("Train Set")
ax[1].legend()

plt.savefig("bdt_output.png", dpi=300)

fpr_test, tpr_test, _ = roc_curve(y_test, y_pred_test)
roc_auc_test = auc(fpr_test, tpr_test)

fpr_train, tpr_train, _ = roc_curve(y_train, y_pred_train)
roc_auc_train = auc(fpr_train, tpr_train)

fig, ax = plt.subplots()
ax.plot(fpr_test, tpr_test, label='Test (area = %0.2f)' % roc_auc_test)
ax.plot(fpr_train, tpr_train, label='Train (area = %0.2f)' % roc_auc_train)

ax.plot([0, 1], [0, 1], 'k--')
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('Background Efficiency')
ax.set_ylabel('Signal Efficiency')
ax.legend()

plt.savefig("roc_curve.png", dpi=300)

features = data.fields
features = [k for k in features if k not in ["weight", "labels"]]

fig, ax = plt.subplots(figsize=(8, 6))
xgb.plot_importance(bdt, ax=ax, height=0.5, importance_type="gain", show_values=False)
ax.set_yticklabels(features) 

r.TMVA.Experimental.SaveXGBoost(bdt, "VBS BDT", "BDT_Weights.root", num_inputs=X.shape[1])




