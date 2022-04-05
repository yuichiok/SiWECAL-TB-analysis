/*Created for eventbuilding development

Before committing changes to the eventbuilding procedure, run the eventbuilding
with both the new and the previous code versions on some example data.
If the trees are different, be sure that this is what you wanted to achieve:

root -l -q -b AssertTreeEquality.C\(\"previous.root\",\"new.root\"\)
*/
Int_t *SortedIndex(TTree *tree, TString tree_name) {
  Int_t n_entries = tree->GetEntries();
  Int_t *index = new Int_t[n_entries];
  if (!tree_name.CompareTo("ecal")) {
    tree->Draw("cycle:event", "", "goff");
    TMath::Sort(n_entries, tree->GetV1(), index, false);
  } else {
    throw std::runtime_error(
        TString::Format("Tree sorting not implemented: %s", tree_name.Data()));
  }
  return index;
}

Int_t TreeEntries(TTree *tree1, TTree *tree2, TString tree_name) {
  Int_t n1 = tree1->GetEntries();
  Int_t n2 = tree2->GetEntries();
  // return (n1 < n2) ? n1 : n2; // Uncomment to allow different-length trees.
  if (n1 != n2) {
    throw std::runtime_error(TString::Format(
        "ERROR: Different number of entries in %s tree: %i vs %i",
        tree_name.Data(), n1, n2));
  }
  return n1;
}

Int_t CompareTrees(TTree *tree1, TTree *tree2, Int_t *id_sorted1,
                   Int_t *id_sorted2, TString tree_name = "ecal",
                   Int_t max_print = 0) {
  Int_t n_entries = TreeEntries(tree1, tree2, tree_name);
  TObjArray *leaves1 = tree1->GetListOfLeaves();
  Int_t n_leaves1 = leaves1->GetEntriesFast();
  TObjArray *leaves2 = tree2->GetListOfLeaves();
  Int_t n_leaves2 = leaves2->GetEntriesFast();
  if (n_leaves1 != n_leaves2) {
    throw std::runtime_error(TString::Format(
        "ERROR: Different number of branches in %s tree: %i vs %i",
        tree_name.Data(), n_leaves1, n_leaves2));
  }
  bool is_equal = true;
  Int_t i_different = 0;
  for (Int_t i_entry = 0; i_entry < n_entries; i_entry++) {
    if (max_print != 0) {
      printf("%3d%% (%8d/%d)\r", Int_t(100 * i_entry / n_entries), i_entry,
             n_entries);
    }
    is_equal = true;
    tree1->GetEntry(id_sorted1[i_entry]);
    tree2->GetEntry(id_sorted2[i_entry]);
    for (Int_t i = 0; i < n_leaves1; i++) {
      Int_t n_hits1 = ((TLeaf *)leaves1->UncheckedAt(i))->GetLen();
      Int_t n_hits2 = ((TLeaf *)leaves2->UncheckedAt(i))->GetLen();
      if (n_hits1 != n_hits2) {
        is_equal = false;
        i_different++;
        break;
      }
      for (Int_t i_hit = 0; i_hit < n_hits1; i_hit++) {
        if (((TLeaf *)leaves1->UncheckedAt(i))->GetValue(i_hit) !=
            ((TLeaf *)leaves2->UncheckedAt(i))->GetValue(i_hit)) {
          is_equal = false;
          break;
        }
      }
      if (!is_equal) {
        i_different++;
        break;
      }
    }
    if (!is_equal & (i_different <= max_print)) {
      cout << endl;
      for (Int_t i = 0; i < n_leaves1; i++) {
        cout << "    ";
        cout << ((TLeaf *)leaves1->UncheckedAt(i))->GetFullName() << ": ";
        Int_t n_hits1 = ((TLeaf *)leaves1->UncheckedAt(i))->GetLen();
        Int_t n_hits2 = ((TLeaf *)leaves2->UncheckedAt(i))->GetLen();
        if ((n_hits1 > 1) | (n_hits2 > 1)) {
          cout << "length " << n_hits1 << " <-> " << n_hits2;
          cout << endl << "        ";
          for (Int_t i_hit = 0; i_hit < n_hits1; i_hit++) {
            cout << ((TLeaf *)leaves1->UncheckedAt(i))->GetValue(i_hit) << ", ";
          }
          cout << endl << "        ";
          for (Int_t i_hit = 0; i_hit < n_hits2; i_hit++) {
            cout << ((TLeaf *)leaves2->UncheckedAt(i))->GetValue(i_hit) << ", ";
          }
        } else {
          cout << ((TLeaf *)leaves1->UncheckedAt(i))->GetValue() << " <-> ";
          cout << ((TLeaf *)leaves2->UncheckedAt(i))->GetValue();
        }
        cout << endl;
      }
      cout << endl;
    }
  }
  if (max_print != 0) {
    cout << endl;
  }
  return i_different;
}

Int_t AssertTreeEquality(TString file_name1, TString file_name2,
                         TString tree_name = "ecal", Int_t max_print = 2) {
  TFile f1(file_name1);
  TTree *tree1 = (TTree *)f1.Get(tree_name);
  TFile f2(file_name2);
  TTree *tree2 = (TTree *)f2.Get(tree_name);
  Int_t n_entries = TreeEntries(tree1, tree2, tree_name);
  Int_t *id_sorted1 = SortedIndex(tree1, tree_name);
  Int_t *id_sorted2 = SortedIndex(tree2, tree_name);
  Int_t i_different =
      CompareTrees(tree1, tree2, id_sorted1, id_sorted2, tree_name, max_print);
  return i_different;
}
