import os

import ROOT

def print_from_one_file(filename):
    tree = "Events"

    df = ROOT.RDataFrame(tree, filename)

    # Get all columns
    cols = list(df.GetColumnNames())

    # Display n rows (set n=None or omit to print all, but be careful!)
    display = df.Display(cols, 100)
    display.Print()


filenames = ["onemu_parthiggs_2016H.root", "twomu_parthiggs_2016H.root",
             "onemu_higgsto4mu_2016G.root", "twomu_higgsto4mu_2016G.root",
             "onemu_higgsto4mu_2016H.root", "twomu_higgsto4mu_2016H.root"]

for filename in filenames:

    if os.path.exists(filename):
        print("Printing the Higgs to 4 leptons events from", filename)
        print_from_one_file(filename)
