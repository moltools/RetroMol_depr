#!/usr/bin/env python
import argparse
import joblib

import numpy as np
from rdkit import Chem, RDLogger
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report 

from retromol.chem import mol_to_fingerprint

def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input compounds with classification label.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to model save file.")
    return parser.parse_args()

def main() -> None:
    """
    Driver code.
    """
    RDLogger.DisableLog("rdApp.*")

    args = cli()

    X = []
    y = []

    with open (args.input, "r") as f:
        for line in f:
            smi, cls = line.strip().split(",")
            cls = int(cls)
            
            mol = Chem.MolFromSmiles(smi)
            fp = mol_to_fingerprint(mol, 4, 2048)

            X.append(fp)
            y.append(cls)

    X = np.array(X)
    y = np.array(y)

    print(X.shape, y.shape)

    # Shuffle X and y and split into train and test
    idx = np.arange(len(X))
    np.random.shuffle(idx)
    X = X[idx]
    y = y[idx]
    split = int(0.8 * len(X))

    X_train = X[:split]
    y_train = y[:split]
    X_test = X[split:]
    y_test = y[split:]

    model = RandomForestClassifier(n_estimators=1000, max_depth=20, random_state=0)

    model.fit(X_train, y_train)

    print("Train accuracy:", model.score(X_train, y_train))
    print("Test accuracy:", model.score(X_test, y_test))

    y_pred = model.predict(X_test)
    print(classification_report(y_test, y_pred))

    joblib.dump(model, args.output)

    exit(0)

if __name__ == "__main__":
    main()
