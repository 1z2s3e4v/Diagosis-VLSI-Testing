/**********************************************************************/
/*           Build the fault list;                                    */
/*           Calculate the total fault coverage.                      */
/*           Building the dummy input and output gate list;           */
/*                                                                    */
/*           Author: Bing-Chen (Benson) Wu                            */
/*           update: 01/21/2018                                       */
/*           Revised: Cheng-Yun Hsieh                                 */
/*           modified:02/05/2020   revise fault collapsing            */
/**********************************************************************/

#include "atpg.h"
#include <unordered_map>

string ATPG::f2s(fptr f){
  char buffer[255];
  sprintf(buffer, "%s %s %s SA%d", sort_wlist[f->to_swlist]->name.c_str(), f->node->name.c_str(), (f->io?"GO":"GI"), f->fault_type);
  string s(buffer);
  return s;
}

void ATPG::generate_fault_list() {
  int fault_num;
  wptr w;
  nptr n;
  fptr_s f;
  unordered_map<wptr, int> num_of_eqv_sa0, num_of_eqv_sa1;
  unordered_map<wptr, unordered_set<string> > eqv_sa0s, eqv_sa1s;
  /* walk through every wire in the circuit*/
  for (auto pos : sort_wlist) {
    w = pos;
    n = w->inode.front();

    /* for each gate, create a gate output stuck-at zero (SA0) fault */
    f = move(fptr_s(new(nothrow) FAULT));
    if (f == nullptr) error("No more room!");
    f->tfsp = test_fails;
    f->node = n;

    f->io = GO;     // gate output SA0 fault
    f->fault_type = STUCK0;
    f->to_swlist = w->wlist_index;
    // for AND NOR NOT BUF, their GI fault is equivalent to GO SA0 fault
    // For each fault f, we count how many faults are in the same equivalent class as f.
    //   So that, we can calculate uncollapsed fault coverage.

    f->eqv_fault_num = 1;    // GO SA0 fault itself
    f->eqv_faults.clear(); f->eqv_faults.insert(f2s(f.get()));
    switch (n->type) {
      case AND:
      case BUF:
        for (wptr wptr_ele: w->inode.front()->iwire) {
          f->eqv_fault_num += num_of_eqv_sa0[wptr_ele];   // count how many equivalent faults in gate input
          // f->eqv_faults.insert(eqv_sa0s[wptr_ele].begin(), eqv_sa0s[wptr_ele].end());
          for(string f_str : eqv_sa0s[wptr_ele]){
            char buffer0[255]; sprintf(buffer0, "%s %s %s SA%d", wptr_ele->name.c_str(), wptr_ele->inode.front()->name.c_str(), "GO", STUCK0); string s0(buffer0);
            if(wptr_ele->onode.size() > 1 && f_str == s0){
              char new_buffer0[255]; sprintf(new_buffer0, "%s %s %s SA%d", wptr_ele->name.c_str(), f->node->name.c_str(), "GI", STUCK0); string new_s0(new_buffer0);
              f->eqv_faults.insert(new_s0);
            } else{
              f->eqv_faults.insert(f_str);
            }
          }
        }
        break;
      case NOR:
      case NOT:
        for (wptr wptr_ele: w->inode.front()->iwire) {
          f->eqv_fault_num += num_of_eqv_sa1[wptr_ele];
          // f->eqv_faults.insert(eqv_sa1s[wptr_ele].begin(), eqv_sa1s[wptr_ele].end());
          for(string f_str : eqv_sa1s[wptr_ele]){
            char buffer1[255]; sprintf(buffer1, "%s %s %s SA%d", wptr_ele->name.c_str(), wptr_ele->inode.front()->name.c_str(), "GO", STUCK1); string s1(buffer1);
            if(wptr_ele->onode.size() > 1 && f_str == s1){
              char new_buffer1[255]; sprintf(new_buffer1, "%s %s %s SA%d", wptr_ele->name.c_str(), f->node->name.c_str(), "GI", STUCK1); string new_s1(new_buffer1);
              f->eqv_faults.insert(new_s1);
            } else{
              f->eqv_faults.insert(f_str);
            }
          }
        }
        break;
      case INPUT:
      case OR:
      case NAND:
      case EQV:
      case XOR:
        break;
    }
    //  these are representative fault in the equivalent fault class
    //  For a wire w,  w.onode is a vector of w's fanout nodes.
    //   w->onode.size is the number of fanout nodes.  if w->onode.size >1, w is a fanout stem.
    //   w->onode.front()->type is the gate type at the output of w (assume w has only 1 fanout)
    //   w->onode.front()->type = OR means w is an input of OR gate.
    if (w->onode.size() > 1 || w->onode.front()->type == EQV || w->onode.front()->type == OR
        || w->onode.front()->type == NOR || w->onode.front()->type == XOR || w->onode.front()->type == OUTPUT) {
      num_of_gate_fault += f->eqv_fault_num; // accumulate total uncollapsed faults
      flist_undetect.push_front(f.get()); // initial undetected fault list contains all faults
      flist.push_front(move(f));  // push into the fault list
    } else {
      num_of_eqv_sa0[w] = f->eqv_fault_num;
      eqv_sa0s[w] = f->eqv_faults;
    }

    /* for each gate, create a gate output stuck-at one (SA1) fault */
    f = move(fptr_s(new(nothrow) FAULT));

    if (f == nullptr) error("No more room!");
    f->tfsp = test_fails;
    f->node = n;
    f->io = GO;
    f->fault_type = STUCK1;
    f->to_swlist = w->wlist_index;
    f->eqv_fault_num = 1;
    f->eqv_faults.clear(); f->eqv_faults.insert(f2s(f.get()));
    /* for OR NAND NOT BUF, their GI fault is equivalent to GO SA1 fault */
    switch (n->type) {
      case OR:
      case BUF:
        for (wptr wptr_ele: w->inode.front()->iwire) {
          if (num_of_eqv_sa1.find(wptr_ele) == num_of_eqv_sa1.end())
            cerr << wptr_ele << " is not in hashmap." << endl;
          f->eqv_fault_num += num_of_eqv_sa1[wptr_ele];
          // f->eqv_faults.insert(eqv_sa1s[wptr_ele].begin(), eqv_sa1s[wptr_ele].end());
          for(string f_str : eqv_sa1s[wptr_ele]){
            char buffer1[255]; sprintf(buffer1, "%s %s %s SA%d", wptr_ele->name.c_str(), wptr_ele->inode.front()->name.c_str(), "GO", STUCK1); string s1(buffer1);
            if(wptr_ele->onode.size() > 1 && f_str == s1){
              char new_buffer1[255]; sprintf(new_buffer1, "%s %s %s SA%d", wptr_ele->name.c_str(), f->node->name.c_str(), "GI", STUCK1); string new_s1(new_buffer1);
              f->eqv_faults.insert(new_s1);
            } else{
              f->eqv_faults.insert(f_str);
            }
          }
        }
        break;
      case NAND:
      case NOT:
        for (wptr wptr_ele: w->inode.front()->iwire) {
          if (num_of_eqv_sa0.find(wptr_ele) == num_of_eqv_sa0.end())
            cerr << wptr_ele << " is not in hashmap." << endl;
          f->eqv_fault_num += num_of_eqv_sa0[wptr_ele];
          // f->eqv_faults.insert(eqv_sa0s[wptr_ele].begin(), eqv_sa0s[wptr_ele].end());
          for(string f_str : eqv_sa0s[wptr_ele]){
            char buffer0[255]; sprintf(buffer0, "%s %s %s SA%d", wptr_ele->name.c_str(), wptr_ele->inode.front()->name.c_str(), "GO", STUCK0); string s0(buffer0);
            if(wptr_ele->onode.size() > 1 && f_str == s0){
              char new_buffer0[255]; sprintf(new_buffer0, "%s %s %s SA%d", wptr_ele->name.c_str(), f->node->name.c_str(), "GI", STUCK0); string new_s0(new_buffer0);
              f->eqv_faults.insert(new_s0);
            } else{
              f->eqv_faults.insert(f_str);
            }
          }
        }
        break;
      case INPUT:
      case AND:
      case NOR:
      case EQV:
      case XOR:
        break;
    }
    if (w->onode.size() > 1 || w->onode.front()->type == EQV || w->onode.front()->type == AND
        || w->onode.front()->type == NAND || w->onode.front()->type == XOR || w->onode.front()->type == OUTPUT) {
      num_of_gate_fault += f->eqv_fault_num; // accumulate total fault count
      flist_undetect.push_front(f.get()); // initial undetected fault list contains all faults
      flist.push_front(move(f));  // push into the fault list
    } else {
      num_of_eqv_sa1[w] = f->eqv_fault_num;
      eqv_sa1s[w] = f->eqv_faults;
    }

    /*if w has multiple fanout branches */
    if (w->onode.size() > 1) {
      num_of_eqv_sa0[w] = num_of_eqv_sa1[w] = 1;
      char buffer0[255]; sprintf(buffer0, "%s %s %s SA%d", w->name.c_str(), w->inode.front()->name.c_str(), "GO", STUCK0); string s0(buffer0); 
      eqv_sa0s[w].insert(s0);
      char buffer1[255]; sprintf(buffer1, "%s %s %s SA%d", w->name.c_str(), w->inode.front()->name.c_str(), "GO", STUCK1); string s1(buffer1); 
      eqv_sa1s[w].insert(s1);
      for (nptr nptr_ele: w->onode) {
        /* create SA0 for OR NOR EQV XOR gate inputs  */
        switch (nptr_ele->type) {
          case OUTPUT:
          case OR:
          case NOR:
          case EQV:
          case XOR:
            f = move(fptr_s(new(nothrow) FAULT));
            if (f == nullptr) error("No more room!");
            f->tfsp = test_fails;
            f->node = nptr_ele;
            f->io = GI;
            f->fault_type = STUCK0;
            f->to_swlist = w->wlist_index;
            f->eqv_fault_num = 1;
            f->eqv_faults.clear(); f->eqv_faults.insert(f2s(f.get()));
            /* f->index is the index number of gate input,
               which GI fault is associated with*/
            for (int k = 0; k < nptr_ele->iwire.size(); k++) {
              if (nptr_ele->iwire[k] == w) f->index = k;
            }
            num_of_gate_fault++;
            flist_undetect.push_front(f.get());
            flist.push_front(move(f));
            break;
        }
        /* create SA1 for AND NAND EQV XOR gate inputs  */
        switch (nptr_ele->type) {
          case OUTPUT:
          case AND:
          case NAND:
          case EQV:
          case XOR:
            f = move(fptr_s(new(nothrow) FAULT));
            if (f == nullptr) error("No more room!");
            f->tfsp = test_fails; 
            f->node = nptr_ele;
            f->io = GI;
            f->fault_type = STUCK1;
            f->to_swlist = w->wlist_index;
            f->eqv_fault_num = 1;
            f->eqv_faults.clear(); f->eqv_faults.insert(f2s(f.get()));
            for (int k = 0; k < nptr_ele->iwire.size(); k++) {
              if (nptr_ele->iwire[k] == w) f->index = k;
            }
            num_of_gate_fault++;
            flist_undetect.push_front(f.get());
            flist.push_front(move(f));
            break;
        }
      }
    }
  }
  flist.reverse();
  flist_undetect.reverse();
  /*walk through all faults, assign fault_no one by one  */
  fault_num = 0;
  for (fptr f: flist_undetect) {
    f->eqv_faults.erase(f2s(f)); // remove it self
    f->fault_no = fault_num;
    fault_num++;
    //cout << f->fault_no << f->node->name << ":" << (f->io?"O":"I") << (f->io?9:(f->index)) << "SA" << f->fault_type << endl;
  }

  if (!diag_only)
    fprintf(stdout, "#number of equivalent faults = %d\n", fault_num);
}/* end of generate_fault_list */

/* computing the actual fault coverage */
void ATPG::compute_fault_coverage() {
  double gate_fault_coverage, eqv_gate_fault_coverage;
  int no_of_detect, eqv_no_of_detect, eqv_num_of_gate_fault;
  fptr f;

  debug = 0;
  no_of_detect = 0;
  gate_fault_coverage = 0.0;
  eqv_no_of_detect = 0;
  eqv_gate_fault_coverage = 0.0;
  eqv_num_of_gate_fault = 0;

  for (const auto &pos : flist) {
    f = pos.get();
    if (debug) {
      if (f->detect != 0) {
        switch (f->node->type) {
          case INPUT:
            fprintf(stdout, "primary input: %s\n", f->node->owire.front()->name.c_str());
            break;
          case OUTPUT:
            fprintf(stdout, "primary output: %s\n", f->node->iwire.front()->name.c_str());
            break;
          default:
            fprintf(stdout, "gate: %s ;", f->node->name.c_str());
            if (f->io == GI) {
              fprintf(stdout, "input wire name: %s\n", f->node->iwire[f->index]->name.c_str());
            } else {
              fprintf(stdout, "output wire name: %s\n", f->node->owire.front()->name.c_str());
            }
            break;
        }
        fprintf(stdout, "fault_type = ");
        switch (f->fault_type) {
          case STUCK0:
            fprintf(stdout, "s-a-0\n");
            break;
          case STUCK1:
            fprintf(stdout, "s-a-1\n");
            break;
        }
        fprintf(stdout, "no of equivalent fault = %d\n", f->eqv_fault_num);
        fprintf(stdout, "detection flag = %d\n", f->detect);
        fprintf(stdout, "\n");
      }
    }
    if (f->detect == TRUE) {
      no_of_detect += f->eqv_fault_num;
      eqv_no_of_detect++;
    }
    eqv_num_of_gate_fault++;
  }
  if (num_of_gate_fault != 0)
    gate_fault_coverage = (((double) no_of_detect) / num_of_gate_fault) * 100;
  if (eqv_num_of_gate_fault != 0)
    eqv_gate_fault_coverage = (((double) eqv_no_of_detect) / eqv_num_of_gate_fault) * 100;

  /* print out fault coverage results */
  fprintf(stdout, "\n");
  fprintf(stdout, "#FAULT COVERAGE RESULTS :\n");
  fprintf(stdout, "#number of test vectors = %d\n", in_vector_no);
  fprintf(stdout,
          "#total number of faults (uncollapsed) = %d\n",
          num_of_gate_fault);  // uncollapsed gate-level fault
  fprintf(stdout, "#total number of detected faults (uncollapsed)= %d\n", no_of_detect);
  fprintf(stdout, "#fault coverage (uncollapsed) = %5.2f%%\n", gate_fault_coverage);  // uncollapsed fault coverage
  fprintf(stdout, "#number of equivalent faults (collapsed) = %d\n", eqv_num_of_gate_fault);
  fprintf(stdout, "#number of equivalent detected faults (collapsed) = %d\n", eqv_no_of_detect);
  fprintf(stdout, "#equivalent fault coverage (collapsed) = %5.2f%%\n", eqv_gate_fault_coverage);
  fprintf(stdout, "\n");
}/* end of compute_fault_coverage */

/* for each PI and PO in the whole circuit,
   create a dummy PI gate to feed each PI wire 
   create a dummy PO gate to feed each PO wire. 
   The reason we need dummy gate is to help to deal with boundary conditions easily.
   */

//
// for a primary input, we have this schematic
//   node1 (dummy PI gate) - wireA - node2 (gate)
// we have two fault on wireA, associated with "node1 output faults" 
//
// for a primary output, we have this schematic
//   node1 (gate) - wireA - node2 (dummy PO gate) 
// we have two fault on wireA, associated with "node1 output faults" 
//

void ATPG::create_dummy_gate() {
  int i;
  int num_of_dummy;
  nptr n;
  char sgate[25];
  char intstring[25];

  num_of_dummy = 0;

  /* create a dummy PI gate for each PI wire */
  for (i = 0; i < cktin.size(); i++) {
    num_of_dummy++;

    /* the dummy gate name is  "dummay_gate#"  */
    sprintf(intstring, "%d", num_of_dummy);
    sprintf(sgate, "dummy_gate%s", intstring);

    /* n is the dummy PI gate, cktin[i]is the original PI wire.  feed n to cktin[i] */
    n = getnode(string(sgate));
    n->type = INPUT;
    n->owire.push_back(cktin[i]);
    cktin[i]->inode.push_back(n);
  } // for i

  /* create a dummy PO gate for each PO wire */
  for (i = 0; i < cktout.size(); i++) {
    num_of_dummy++;

    /* the dummy gate name is  "dummay_gate#"  */
    sprintf(intstring, "%d", num_of_dummy);
    sprintf(sgate, "dummy_gate%s", intstring);

    /* n is the dummy PO gate, cktout[i] is the original PO wire.  feed cktout[i] to n */
    n = getnode(sgate);
    n->type = OUTPUT;
    n->iwire.push_back(cktout[i]);
    cktout[i]->onode.push_back(n);
  } // for i
}/* end of create_dummy_gate */

char ATPG::itoc(const int &i) {
  switch (i) {
    case 1:
      return '1';
    case 2:
      return 'U';
    case 0:
      return '0';
  }
}
