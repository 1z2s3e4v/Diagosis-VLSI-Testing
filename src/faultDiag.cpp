/**********************************************************************/
/*  Parallel-Fault Event-Driven Transition Delay Fault Simulator      */
/*                                                                    */
/*           Author: Tsai-Chieh Chen                                  */
/*           last update : 10/22/2018                                 */
/**********************************************************************/

#include "atpg.h"
#include <unordered_set>
#define num_of_faults_in_parallel  16

/* The faulty_wire contains a list of wires that 
 * change their values in the fault simulation for a particular packet.
 * (that is wire_value_g != wire_value_f) 
 * Note that the wire themselves are not necessarily a fault site.
 * The list is linked by the pnext pointers */
string failLog_str[2] = {"expect L, observe H", "expect H, observe L"};

#include <bitset>
#include <string>
std::string toBinary(uint32_t num) {
    std::bitset<32> binary(num);
    return binary.to_string();
}

void ATPG::print_circuit_summary(){
    fprintf(stdout, "\n");
    fprintf(stdout, "#Circuit Summary:\n");
    fprintf(stdout, "#---------------\n");
    fprintf(stdout, "#number of inputs = %d\n", int(cktin.size()));
    fprintf(stdout, "#number of outputs = %d\n", int(cktout.size()));
    fprintf(stdout, "#number of gates = %d\n", ncktnode);
    fprintf(stdout, "#number of wires = %d\n", ncktwire);
    fprintf(stdout, "#number of vectors = %d\n", tr_unexamined1.size()); 
    fprintf(stdout, "#number of failing outputs = %d\n", test_fails);
    fprintf(stdout, "\n");
}

//void ATPG::construct_po_map(){ // TODO: Step 0. 
    
//}

void ATPG::ranking(){
    int flist_size = ranks.size();
    int i, j;
    bool swapped;
    for (i = 0; i < flist_size - 1; ++i) {
        swapped = false;
        for (j = 0; j < flist_size - i -1 ; ++j) {
            if (ranks[j]->score < ranks[j + 1]->score) {
                swap(ranks[j], ranks[j + 1]);
                swapped = true;
            }
        }
        if (swapped == false) break;
    }
}

/* fault simulate a set of test vectors */
void ATPG::genFailLog_fault_simulation(int &total_detect_num) {
  int i;
  int current_detect_num = 0;
  
  /* for every vector */
  for (i = vectors.size() - 1; i >= 0; i--) {
    genFailLog_fault_sim_a_vector(vectors[i], current_detect_num, i); // print fail log
    total_detect_num += current_detect_num;
  }
}// genFailLog_fault_simulation
/* diagnosis */
void ATPG::diag(){
  
  string  wire_name;
  int observed, expected;
  trptr fc;
  trptr temp;
  int flist_size;
  for (auto pos = flist_undetect.cbegin(); pos != flist_undetect.cend(); ++pos) {
    flist_undetect1.push_back(*pos);
    SF.push_back(*pos);
    
  }

  
  for (auto vec : vectors){
    if(tr_unexamined1.find(vec) != tr_unexamined1.end())  FT_set.push_back(vec);  //failing test set
    else PT_set.push_back(vec); //passing test set
  }
  
  flist_size = flist_undetect1.size();
  //Phase 1
  for (auto vec : PT_set) {
    single_SAF_simulation(vec);
  }

  for (auto f : SF)
  {
    SF1[sort_wlist[(f)->to_swlist]->name].push_back(f);
  }
  for (auto f : PNF)
  {
    PNF1[sort_wlist[(f)->to_swlist]->name].push_back(f);
  }

  int prev = -1;
  bool phase4 = false;
  for (int timer = 0; timer < 20  ; timer++)
  {
    int i = 0;
    int drop1 =0, drop2 =0;
    //cout << timer << " " << SF.size() << " " << PNF.size()  <<endl;
    //Phase 2
    for(int fi = 0 ; fi < FT_set.size(); fi++)
    {
      //cout << "\n phase 2 " <<FT_set[fi] << endl;
      for( auto f : SF)
      {
        //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;
      }
      multiple_SAF_simulation(FT_set[fi]);
      if(!mismatching_output.empty()){
        //for(auto m : mismatching_output) cout << m << endl;
        single_SAF_simulation23(FT_set[fi], drop1);i++;
      
      }
      //if(timer%3 == 0) reverse(SF.begin(), SF.end());
      mismatching_output.erase(mismatching_output.begin(), mismatching_output.end());
      for( auto f : SF)
      {
        //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;
      }
      multiple_SAF_simulation(FT_set[fi]);
      if(!mismatching_output.empty()){
        //for(auto m : mismatching_output) cout << m << endl;
        single_SAF_simulation23(FT_set[fi], drop1);i++;
      
      }
      mismatching_output.erase(mismatching_output.begin(), mismatching_output.end());
    }
    //Phase 3
    for(int pi = 0 ; pi < PT_set.size(); pi++)
    {
      //cout << "\n phase 3 " <<PT_set[pi]  << endl;
      for( auto f : SF)
      {
        //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;
      }
      multiple_SAF_simulation1(PT_set[pi]);
      //cout << timer << " " << SF.size() << " " << PNF.size()  <<endl;
      if(!mismatching_output.empty()){
        //for(auto m : mismatching_output) cout << m << endl;
        single_SAF_simulation34(PT_set[pi], drop2);i++ ;
      }
      //if(timer%5 == 0) reverse(SF.begin(), SF.end());
      mismatching_output.erase(mismatching_output.begin(), mismatching_output.end());
      for( auto f : SF)
      {
        //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;
      }
      multiple_SAF_simulation1(PT_set[pi]);
      //cout << timer << " " << SF.size() << " " << PNF.size()  <<endl;
      if(!mismatching_output.empty()){
        //for(auto m : mismatching_output) cout << m << endl;
        single_SAF_simulation34(PT_set[pi], drop2);i++ ;
      }
      mismatching_output.erase(mismatching_output.begin(), mismatching_output.end());
    }
    reverse(SF.begin(), SF.end());
    
    if ((prev == SF.size()) &&  (i == 0) && (timer%2 == 0) ){//cout << "con1 con2 "  << timer << endl; 
    phase4 = true;break;}
    prev = SF.size();
    //cout << SF.size() <<  endl;
  }
  for (int i = 0; i < sort_wlist.size(); i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed = FALSE;
  }
  for( auto f : SF){
    f->detect = FALSE;
  }
  //cout << endl;
  diag1();
  
  std::sort(SF.begin(), SF.end(),[&](const fptr fptr_ele1, const fptr fptr_ele2) {return fptr_ele1->score < fptr_ele2->score;} );
  for( auto f : SF)
  {
    //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << " " << f->score << endl;
  }
  //cout << endl;
  
  // Phase 4
  bool flag = false;
  int cc = SF.size();
  if(phase4){
    while(cc--){
      bool flag1 = false;
      for(auto f : SF)
      {
        bool flag2 = true;
        if(flag1 == true) continue;
        for (auto vec : vectors){
          if (flag2 == false) continue;
          flag2 = multiple_SAF_simulation2(vec, f);
          //cout << vec << " " << flag2 <<" "<<f->fault_no  <<endl;
        }
        if(flag2==true) // fault is removable
        {
          f->remove = true;
          flag1 = true;
        }
      }
    
    SF.erase(
              std::remove_if(SF.begin(), SF.end(),
                  [&](const fptr fptr_ele) {
                  if ((fptr_ele->remove == TRUE )&& (fptr_ele->score != 100)) {
                    //cout <<"drop6: "  <<fptr_ele->fault_no << " " << fptr_ele->node->name << ":" << (fptr_ele->io?"O":"I")<< " "  << sort_wlist[fptr_ele->to_swlist]->name << "SA" << fptr_ele->fault_type << endl;
                    //cout << fptr_ele->fault_no  << " " <<fptr_ele->score <<  endl;
                    return true;
                  } else if (fptr_ele->remove == FALSE){
                    return false;
                  }
                }),
              SF.end());
    }
  }
  for(auto f : SF)
  {
    //f->score = ((double) f->tfsf*10)/((double) (f->tfsf*10 + f->tfsp*10 + f->tpsf))*100;
    if(f->score > 0)  ranks.push_back(f);
  }
  //cout << SF.size() << endl;
  for( auto f : SF)
  {
    //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;
  }
  //diag1();
}
/* fault simulate a single test vector */
void ATPG::genFailLog_fault_sim_a_vector(const string &vec, int &num_of_current_detect, int &vec_index) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;
  set<wptr> failed_ws;

  num_of_fault = 0; // counts the number of faults in a packet

  /* num_of_current_detect is used to keep track of the number of undetected faults
   * detected by this vector.  Initialize it to zero */
  num_of_current_detect = 0;

  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;
 
  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) */
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); /* do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i
  /* walk through every undetected fault
   * the undetected fault list is linked by pnext_undetect */
  for (auto pos = examined_faults.cbegin(); pos != examined_faults.cend(); ++pos) {
    f = *pos;
    if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */

    /* if f is a primary output or is directly connected to an primary output
      * the fault is detected */
    if ((f->node->type == OUTPUT) ||
        (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
      /* if this wire is not yet marked as faulty, mark the wire as faulty
        * and insert the corresponding wire to the list of faulty wires. */
      if (!(sort_wlist[f->to_swlist]->is_faulty())) {
        sort_wlist[f->to_swlist]->set_faulty();
        wlist_faulty.push_front(sort_wlist[f->to_swlist]);
      }
      inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
      /* mark the wire as having a fault injected
        * and schedule the outputs of this gate */
      sort_wlist[f->to_swlist]->set_fault_injected();
      /* increment the number of simulated faults in this packet */
      /* start_wire_index keeps track of the smallest level of fault in this packet.
        * this saves simulation time.  */
      start_wire_index = min(start_wire_index, f->to_swlist);
    } else {
      /* if f is an gate output fault */
      if (f->io == GO) {
        /* if this wire is not yet marked as faulty, mark the wire as faulty
          * and insert the corresponding wire to the list of faulty wires. */
        if (!(sort_wlist[f->to_swlist]->is_faulty())) {
          sort_wlist[f->to_swlist]->set_faulty();
          wlist_faulty.push_front(sort_wlist[f->to_swlist]);
        }

        /* add the fault to the simulated fault list and inject the fault */
        inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);

        /* mark the wire as having a fault injected
          * and schedule the outputs of this gate */
        sort_wlist[f->to_swlist]->set_fault_injected();
        for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
          pos_n->owire.front()->set_scheduled();
        }

        /* increment the number of simulated faults in this packet */
        /* start_wire_index keeps track of the smallest level of fault in this packet.
          * this saves simulation time.  */
        start_wire_index = min(start_wire_index, f->to_swlist);
      }  // if gate output fault

      /* the fault is a gate input fault */
      else {
        /* if the fault is propagated, set faulty_wire equal to the faulty wire.
          * faulty_wire is the gate output of f.  */
        faulty_wire = get_faulty_wire(f, fault_type);
        if (faulty_wire != nullptr) {
          /* if the faulty_wire is a primary output, it is detected */
          if (faulty_wire->is_output()) {
            f->detected_time++; /// BONUS
            if(f->detected_time >= detected_num){
              f->detect = TRUE;
            }
          } else {
            /* if faulty_wire is not already marked as faulty, mark it as faulty
              * and add the wire to the list of faulty wires. */
            if (!(faulty_wire->is_faulty())) {
              faulty_wire->set_faulty();
              wlist_faulty.push_front(faulty_wire);
            }
            if(!sort_wlist[f->to_swlist]->is_fault_injected())
            {
              /* add the fault to the simulated list and inject it */
              inject_fault_value(faulty_wire, num_of_fault, fault_type);
              /* mark the faulty_wire as having a fault injected
              *  and schedule the outputs of this gate */
              faulty_wire->set_fault_injected();
            }
            for (auto pos_n : faulty_wire->onode) {
              pos_n->owire.front()->set_scheduled();
            }

            start_wire_index = min(start_wire_index, f->to_swlist);
          }
        }
      }
    } // if  gate input fault

    
  } // end loop. for f = flist

  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f << endl;
      fault_sim_evaluate(sort_wlist[i]);
    }
  } /* event evaluations end here */
  
  while (!wlist_faulty.empty()) {
    w = wlist_faulty.front();
    wlist_faulty.pop_front();
    w->remove_faulty();
    w->remove_fault_injected();
    w->set_fault_free();

    if(w->is_output()){
      // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
      int detected = w->wire_value_g ^ w->wire_value_f;
      for(int f_idx=0; f_idx<1; ++f_idx){ // for every undetected fault
        // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
        // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
        if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
          failed_ws.insert(w);
          
        }
      }
      //printf("%s %x \n",w->name.c_str(),  w->wire_value_f);
    }
  } // pop out all faulty wires
  num_of_fault = 0;  // reset the counter of faults in a packet
  start_wire_index = 10000;  //reset this index to a very large value.
  for(wptr failed_w : failed_ws){
    //printf("vector[%d] %s %s  # T'%s'\n", vec_index, failed_w->name.c_str(), failLog_str[failed_w->wire_value_g].c_str(), vec.c_str());
    //cout << failed_w->wire_value_g << endl;
    if(failed_w->value == 1){
      printf("vector[%d] %s %s  # T'%s' \n", vec_index, failed_w->name.c_str(), failLog_str[1].c_str(), vec.c_str());
    } else if (failed_w->value == 0){
      printf("vector[%d] %s %s  # T'%s' \n", vec_index, failed_w->name.c_str(), failLog_str[0].c_str(), vec.c_str());
    }
  }
}/* end of fault_sim_a_vector */

void ATPG::set_examined_faults(const string &wire, const string &gate, const string &io, const string &fault_type) {
//   cout << wire << " " << gate << " " << io << " " <<fault_type << endl;
  int fault_num;
  wptr w;
  nptr n;
  fptr_s f;
  bool flag= 0;
  unordered_map<wptr, int> num_of_eqv_sa0, num_of_eqv_sa1;
  for (wptr pos : sort_wlist) {
    if((pos->name != wire ) || ((io == "GO") && (pos->inode.front()->name != gate))) continue;
    if(io == "GI") {
      for(nptr nn : pos->onode) {  
        if (nn->name == gate) {
          flag = 1;
          n = nn;
        }
      }
      if (!flag) continue;
    }
    else if(io == "GO") {
      n = pos->inode.front();
    }
    w = pos;
    f = move(fptr_s(new (nothrow) FAULT));
    if (f == nullptr)
      error("No more room!");
    f->node = n;
    //cout << "dsad " << n->name << endl;
    f->io = io == "GO"; // gate output SA0 fault
    f->fault_type = fault_type == "SA1" ? STUCK1:STUCK0 ;
    f->to_swlist = w->wlist_index;
    examined_faults.push_front(f.get());     // initial undetected fault list contains all faults
    flist.push_front(move(f));             // push into the fault list
  }
  
  //flist.reverse();
  examined_faults.reverse();
  /*walk through all faults, assign fault_no one by one  */
  fault_num = 0;
  for (fptr f : examined_faults) {
    f->fault_no = fault_num;
    fault_num++;
    //cout << f->fault_no << " "  << f->node->name << ":" << (f->io?"O":"I") <<" "  << (f->io?9:(f->index)) <<" "  << "SA" << f->fault_type << endl;
  }
  // fprintf(stdout, "#number of equivalent faults = %d\n", fault_num);
}
void ATPG::MSAF_simulation(vector<int> &bach)
{

}

/* fault simulate a single test vector */
void ATPG::single_SAF_simulation(const string &vec) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) */
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); /* do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i

  /* walk through every undetected fault
   * the undetected fault list is linked by pnext_undetect */
  for (int index = 0; index < SF.size(); index ++) {
    f = SF[index];
    f->detect = FALSE;
    if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */

    /* consider only active (aka. excited) fault
     * (sa1 with correct output of 0 or sa0 with correct output of 1) */
    if (f->fault_type != sort_wlist[f->to_swlist]->value) {

      /* if f is a primary output or is directly connected to an primary output
       * the fault is detected */
      if ((f->node->type == OUTPUT) ||
          (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
        f->detected_time++; /// BONUS
        if(f->detected_time >= detected_num){
          f->detect = TRUE;
        }
        if(!tr_unexamined1[vec].empty())
        {
          if((tr_unexamined1[vec].find(sort_wlist[f->to_swlist]->name) != tr_unexamined1[vec].end())) { 
            f->tfsf++;
            f->tfsp--;
          }
          else {f->tpsf++;}
          // no fault is supposed to occur
        } else {f->tpsf++;}
      } else {

        /* if f is an gate output fault */
        if (f->io == GO) {

          /* if this wire is not yet marked as faulty, mark the wire as faulty
           * and insert the corresponding wire to the list of faulty wires. */
          if (!(sort_wlist[f->to_swlist]->is_faulty())) {
            sort_wlist[f->to_swlist]->set_faulty();
            wlist_faulty.push_front(sort_wlist[f->to_swlist]);
          }

          /* add the fault to the simulated fault list and inject the fault */
          simulated_fault_list[num_of_fault] = f;
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);

          /* mark the wire as having a fault injected
           * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
            pos_n->owire.front()->set_scheduled();
          }

          /* increment the number of simulated faults in this packet */
          num_of_fault++;
          /* start_wire_index keeps track of the smallest level of fault in this packet.
           * this saves simulation time.  */
          start_wire_index = min(start_wire_index, f->to_swlist);
        }  // if gate output fault

          /* the fault is a gate input fault */
        else {

          /* if the fault is propagated, set faulty_wire equal to the faulty wire.
           * faulty_wire is the gate output of f.  */
          faulty_wire = get_faulty_wire(f, fault_type);
          if (faulty_wire != nullptr) {

            /* if the faulty_wire is a primary output, it is detected */
            if (faulty_wire->is_output()) {
              f->detected_time++; /// BONUS
              if(f->detected_time >= detected_num){
                f->detect = TRUE;
              }
              if(!tr_unexamined1[vec].empty())
              {
                if((tr_unexamined1[vec].find(faulty_wire->name) != tr_unexamined1[vec].end())) { 
                  f->tfsf++;
                  f->tfsp--;
                }
                else {f->tpsf++;}
                // no fault is supposed to occur
              } else {f->tpsf++;}
            } else {
              /* if faulty_wire is not already marked as faulty, mark it as faulty
               * and add the wire to the list of faulty wires. */
              if (!(faulty_wire->is_faulty())) {
                faulty_wire->set_faulty();
                wlist_faulty.push_front(faulty_wire);
              }

              /* add the fault to the simulated list and inject it */
              simulated_fault_list[num_of_fault] = f;
              inject_fault_value(faulty_wire, num_of_fault, fault_type);

              /* mark the faulty_wire as having a fault injected
               *  and schedule the outputs of this gate */
              faulty_wire->set_fault_injected();
              for (auto pos_n : faulty_wire->onode) {
                pos_n->owire.front()->set_scheduled();
              }

              num_of_fault++;
              start_wire_index = min(start_wire_index, f->to_swlist);
            }
          }
        }
      } // if  gate input fault
    } // if fault is active


    /*
     * fault simulation of a packet
     */

    /* if this packet is full (16 faults)
     * or there is no more undetected faults remaining (pos points to the final element of flist_undetect),
     * do the fault simulation */
    if ((num_of_fault == num_of_faults_in_parallel) || (index == SF.size() - 1)) {

      /* starting with start_wire_index, evaulate all scheduled wires
       * start_wire_index helps to save time. */
      for (i = start_wire_index; i < nckt; i++) {
        if (sort_wlist[i]->is_scheduled()) {
          sort_wlist[i]->remove_scheduled();
          fault_sim_evaluate(sort_wlist[i]);
        }
      } /* event evaluations end here */

      /* pop out all faulty wires from the wlist_faulty
        * if PO's value is different from good PO's value, and it is not unknown
        * then the fault is detected.
        *
        * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
      */
      while (!wlist_faulty.empty()) {
        w = wlist_faulty.front();
        wlist_faulty.pop_front();
        w->remove_faulty();
        w->remove_fault_injected();
        w->set_fault_free();
        /*TODO*/
        /*
         * After simulation is done,if wire is_output(), we should compare good value(wire_value_g) and faulty value(wire_value_f). 
         * If these two values are different and they are not unknown, then the fault is detected.  We should update the simulated_fault_list.  Set detect to true if they are different.
         * Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
         * After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g).
         */
        // 1.1. if wire is_output()
        if(w->is_output()){
          // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
          int detected = w->wire_value_g ^ w->wire_value_f;
          for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for each faults
            // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
            // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
            if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              simulated_fault_list[f_idx]->detect = TRUE;
              if(!tr_unexamined1[vec].empty())
                {
                  if((tr_unexamined1[vec].find(w->name) != tr_unexamined1[vec].end())) { 
                    simulated_fault_list[f_idx]->tfsf++;
                    simulated_fault_list[f_idx]->tfsp--;
                  }
                  else {simulated_fault_list[f_idx]->tpsf++;}
                  // no fault is supposed to occur
                } else {simulated_fault_list[f_idx]->tpsf++;}
            }
          }
        }
        // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
        w->wire_value_f = w->wire_value_g;
        /*end of TODO*/
        for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ /// BONUS
          if(simulated_fault_list[f_idx]->detect == TRUE){
            simulated_fault_list[f_idx]->detected_time++;
            if(simulated_fault_list[f_idx]->detected_time < detected_num){
              simulated_fault_list[f_idx]->detect = FALSE;
            }
          }
        }

      } // pop out all faulty wires
      num_of_fault = 0;  // reset the counter of faults in a packet
      start_wire_index = 10000;  //reset this index to a very large value.
    } // end fault sim of a packet
  } // end loop. for f = flist
  SF.erase(
    std::remove_if(SF.begin(), SF.end(),
        [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          //fptr_ele->detect == FALSE;
          PNF.push_back(fptr_ele);
          PNF1[sort_wlist[fptr_ele->to_swlist]->name].push_back(fptr_ele);
          return true;
        } else {
          return false;
        }
      }),
    SF.end());
  /* fault dropping  */

}/* end of fault_sim_a_vector */

/* multiple SAF simulate a single test vector */
void ATPG::multiple_SAF_simulation(const string &vec) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) */
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); /* do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed = FALSE;
    
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    f->detect = false;
    if ((f->node->type == OUTPUT) ||
        (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
      if (!(sort_wlist[f->to_swlist]->is_faulty())) {
        sort_wlist[f->to_swlist]->set_faulty();
        wlist_faulty.push_front(sort_wlist[f->to_swlist]);
      }
      if (sort_wlist[f->to_swlist]->fixed == FALSE){ 
        inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
        /* mark the wire as having a fault injected
          * and schedule the outputs of this gate */
        sort_wlist[f->to_swlist]->set_fault_injected();
        sort_wlist[f->to_swlist]->fixed = TRUE;
        start_wire_index = min(start_wire_index, f->to_swlist);
      } 
      else
        {
          //f->detect =TRUE;
        }
    }
    else 
    {
      if (f->io == GO) {
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f << " " << f->fault_type  << " fixed? " << sort_wlist[f->to_swlist]->fixed<<endl;
        if (!(sort_wlist[f->to_swlist]->is_faulty())) {
          sort_wlist[f->to_swlist]->set_faulty();
          wlist_faulty.push_front(sort_wlist[f->to_swlist]);
        }
        if (sort_wlist[f->to_swlist]->fixed == FALSE){
          //cout << " inject "<< endl; 
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
          /* mark the wire as having a fault injected
            * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          sort_wlist[f->to_swlist]->fixed = TRUE;
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
        else
        {
          //f->detect =TRUE;
        }
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f <<endl;
      }
      else
      { 
        continue;
      }
    }
  }
  for (i = 0; i < 0; i++) {
    printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
  }
  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  } /* event evaluations end here */

  
  start_wire_index = 10000;
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    if ((f->node->type == OUTPUT) ||
        (f->io == GO )) {
        continue;
    }
    else // Gate input
    {
      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f1 = wire->wire_value_f;
            inject_fault_value(wire, num_of_fault, ff->fault_type);
            
          }
        }
      }
      //cout << sort_wlist[ f->to_swlist]->name.c_str() << sort_wlist[ f->to_swlist]->wire_value_f1 << " " << sort_wlist[ f->to_swlist]->wire_value_f << endl;

      sort_wlist[ f->to_swlist]->wire_value_f1 = sort_wlist[ f->to_swlist]->wire_value_f;

      faulty_wire = get_faulty_wire(f, fault_type);
      //faulty_wire = get_faulty_wire(f, fault_type);
      if (faulty_wire != nullptr){
        sort_wlist[ f->to_swlist]->wire_value_f1 = sort_wlist[ f->to_swlist]->wire_value_f;
        inject_fault_value( sort_wlist[ f->to_swlist], num_of_fault, f->fault_type);
        faulty_wire = f->node->owire.front();
        
        if(faulty_wire->fixed == FALSE) fault_sim_evaluate(faulty_wire);
        else {
          //cout <<  sort_wlist[ f->to_swlist]->name << " " << faulty_wire->name << faulty_wire->fixed << endl;
          //f->detect = TRUE;
        }
        sort_wlist[ f->to_swlist]->wire_value_f = sort_wlist[ f->to_swlist]->wire_value_f1;
      }
      else
      {
        //cout <<  sort_wlist[ f->to_swlist]->name << " null " << endl;
      }
      
      sort_wlist[ f->to_swlist]->wire_value_f = sort_wlist[ f->to_swlist]->wire_value_f1;
      //cout << sort_wlist[ f->to_swlist]->name.c_str() << sort_wlist[ f->to_swlist]->wire_value_f1 << " " << sort_wlist[ f->to_swlist]->wire_value_f << endl;

      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f = wire->wire_value_f1;
          }
        }
      }
      if (faulty_wire != nullptr ) {
        /* if the faulty_wire is a primary output, it is detected */
        if (faulty_wire->fixed == FALSE ) {
          if (!(faulty_wire->is_faulty())) {
            faulty_wire->set_faulty();
            wlist_faulty.push_front(faulty_wire);
          }
          /* add the fault to the simulated list and inject it */
          //inject_fault_value(faulty_wire, num_of_fault, fault_type);
          /* mark the faulty_wire as having a fault injected
          *  and schedule the outputs of this gate */
          faulty_wire->set_fault_injected();
          faulty_wire->fixed == TRUE;
          if(!faulty_wire->is_output()){
            for (auto pos_n : faulty_wire->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
            }
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
      }
    }
  }

  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  }

  /* pop out all faulty wires from the wlist_faulty
    * if PO's value is different from good PO's value, and it is not unknown
    * then the fault is detected.
    *
    * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
  */
  for (i = 0; i < 0; i++) {
    printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
  }
  for (auto ww : cktout) {
    //printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
    w = ww;
    if(w->is_output()){
      // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
      int detected = w->wire_value_g ^ w->wire_value_f;
      //cout << w->name << vec << " "<< detected << endl;
      for(int f_idx=0; f_idx<1; ++f_idx){ // for each faults
        // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
        // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
        if(((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
          // If sim result is different from observed value -> record mismatching output
          //cout <<w->name << w->wire_value_f << w->wire_value_g <<endl;
          if(((Mask[f_idx] & detected) != 0)&&(tr_unexamined1[vec].find(w->name) == tr_unexamined1[vec].end())) {
            //cout << w->name << vec << " "<< detected << endl;
            mismatching_output.insert(w->name);
          }
          else if (((Mask[f_idx] & detected) == 0)&&(tr_unexamined1[vec].find(w->name) != tr_unexamined1[vec].end())){
            //cout << vec<< w->name  << endl;
            mismatching_output.insert(w->name);
          }
        }

      }
    }
  }
  while (!wlist_faulty.empty()) {
    w = wlist_faulty.front();
    wlist_faulty.pop_front();
    w->remove_faulty();
    w->remove_fault_injected();
    w->set_fault_free();
    // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
    w->wire_value_f = w->wire_value_g;


  } // pop out all faulty wires

  //for (auto m :mismatching_output ) cout << m << endl;
  num_of_fault = 0;  // reset the counter of faults in a packet
  start_wire_index = 10000;  //reset this index to a very large value.
  SF.erase(
    std::remove_if(SF.begin(), SF.end(),
        [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          //cout <<"drop6: "  <<fptr_ele->fault_no << " " << fptr_ele->node->name << ":" << (fptr_ele->io?"O":"I")<< " "  << sort_wlist[fptr_ele->to_swlist]->name << "SA" << fptr_ele->fault_type << endl;
          PNF.push_back(fptr_ele);
          //cout << sort_wlist[fptr_ele->to_swlist]->name << " " <<fptr_ele->detect <<  endl;

          return true;
        } else if (fptr_ele->detect == FALSE){
          return false;
        }
      }),
    SF.end());
  
  for(auto v : SF1)
  {
    for (auto it = SF1[v.first].begin(); it != SF1[v.first].end();)
    {
      fptr fptr_ele1 = *it;
      if (fptr_ele1->detect == TRUE) {
        //cout <<"drop11: "  <<fptr_ele1->fault_no << " " << fptr_ele1->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[fptr_ele1->to_swlist]->name << "SA" << fptr_ele1->fault_type << endl;

        if(std::find(PNF1[v.first].begin(), PNF1[v.first].end(), fptr_ele1) == PNF1[v.first].end() ) 
        {
          PNF1[v.first].push_back(fptr_ele1);
        }
        it = SF1[v.first].erase(it);
      }
      else ++it;
    }
  }
}/* end of fault_sim_a_vector */
/* multiple SAF simulate a single test vector */
void ATPG::multiple_SAF_simulation1(const string &vec) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) */
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); /* do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed = FALSE;
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    f->detect = false;
    if ((f->node->type == OUTPUT) ||
        (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
      if (!(sort_wlist[f->to_swlist]->is_faulty())) {
        sort_wlist[f->to_swlist]->set_faulty();
        wlist_faulty.push_front(sort_wlist[f->to_swlist]);
      }
      if (sort_wlist[f->to_swlist]->fixed == FALSE){ 
        inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
        /* mark the wire as having a fault injected
          * and schedule the outputs of this gate */
        sort_wlist[f->to_swlist]->set_fault_injected();
        sort_wlist[f->to_swlist]->fixed = TRUE;
        start_wire_index = min(start_wire_index, f->to_swlist);
      } else
      {
        f->detect =TRUE;
      }
    }
    else 
    {
      if (f->io == GO) {
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f << " " << f->fault_type  << " fixed? " << sort_wlist[f->to_swlist]->fixed<<endl;
        if (!(sort_wlist[f->to_swlist]->is_faulty())) {
          sort_wlist[f->to_swlist]->set_faulty();
          wlist_faulty.push_front(sort_wlist[f->to_swlist]);
        }
        if (sort_wlist[f->to_swlist]->fixed == FALSE){
          //cout << " inject "<< endl; 
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
          /* mark the wire as having a fault injected
            * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          sort_wlist[f->to_swlist]->fixed = TRUE;
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
        else
        {
          f->detect = true;
        }
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f <<endl;
      }
      else
      { 
        continue;
      }
    }
  }

  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  } /* event evaluations end here */


  start_wire_index = 10000;
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    if ((f->node->type == OUTPUT) ||
        (f->io == GO )) {
        continue;
    }
    else // Gate input
    {
      int tmp = sort_wlist[ f->to_swlist]->wire_value_f;
      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f1 = wire->wire_value_f;
            inject_fault_value(wire, num_of_fault, ff->fault_type);
            
          }
        }
      }
      faulty_wire = get_faulty_wire(f, fault_type);
      if (faulty_wire != nullptr){
        sort_wlist[ f->to_swlist]->wire_value_f1 = sort_wlist[ f->to_swlist]->wire_value_f;
        inject_fault_value( sort_wlist[ f->to_swlist], num_of_fault, f->fault_type);
        faulty_wire = f->node->owire.front();
        //cout <<  sort_wlist[ f->to_swlist]->name << " " << faulty_wire->name << faulty_wire->fixed << endl;
        if(faulty_wire->fixed == FALSE) fault_sim_evaluate(faulty_wire);
        else f->detect = TRUE;
        sort_wlist[ f->to_swlist]->wire_value_f = sort_wlist[ f->to_swlist]->wire_value_f1;
      }
      else
      {
        //cout <<  sort_wlist[ f->to_swlist]->name << " null " << endl;
      }
      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f = wire->wire_value_f1;
          }
        }
      }
      if (faulty_wire != nullptr ) {
        /* if the faulty_wire is a primary output, it is detected */
        if (faulty_wire->fixed == FALSE ) {
          if (!(faulty_wire->is_faulty())) {
            faulty_wire->set_faulty();
            wlist_faulty.push_front(faulty_wire);
          }
          /* add the fault to the simulated list and inject it */
          //inject_fault_value(faulty_wire, num_of_fault, fault_type);
          /* mark the faulty_wire as having a fault injected
          *  and schedule the outputs of this gate */
          faulty_wire->set_fault_injected();
          faulty_wire->fixed == TRUE;
          if(!faulty_wire->is_output()){
            for (auto pos_n : faulty_wire->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
            }
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
      }
    }
  }

  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  }
  for (i = 0; i < 0; i++) {
    printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
  }
  while (!wlist_faulty.empty()) {
    w = wlist_faulty.front();
    wlist_faulty.pop_front();
    w->remove_faulty();
    w->remove_fault_injected();
    w->set_fault_free();
    /*TODO*/
    /*
      * After simulation is done,if wire is_output(), we should compare good value(wire_value_g) and faulty value(wire_value_f). 
      * If these two values are different and they are not unknown, then the fault is detected.  We should update the simulated_fault_list.  Set detect to true if they are different.
      * Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
      * After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g).
      */
    // 1.1. if wire is_output()
    if(w->is_output()){
      // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
      int detected = w->wire_value_g ^ w->wire_value_f;
      for(int f_idx=0; f_idx<1; ++f_idx){ // for each faults
        // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
        // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
        if(((Mask[f_idx] & detected) != 0) &&((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
          mismatching_output.insert(w->name);
        }

      }
    }

    // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
    w->wire_value_f = w->wire_value_g;


  } // pop out all faulty wires

  //for (auto m :mismatching_output ) cout << m << endl;
  num_of_fault = 0;  // reset the counter of faults in a packet
  start_wire_index = 10000;  //reset this index to a very large value.
  SF.erase(
    std::remove_if(SF.begin(), SF.end(),
        [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          //cout <<"drop6: "  <<fptr_ele->fault_no << " " << fptr_ele->node->name << ":" << (fptr_ele->io?"O":"I")<< " "  << sort_wlist[fptr_ele->to_swlist]->name << "SA" << fptr_ele->fault_type << endl;
          PNF.push_back(fptr_ele);
          //cout << sort_wlist[fptr_ele->to_swlist]->name << " " <<fptr_ele->detect <<  endl;

          return true;
        } else if (fptr_ele->detect == FALSE){
          return false;
        }
      }),
    SF.end());
  
  for(auto v : SF1)
  {
    for (auto it = SF1[v.first].begin(); it != SF1[v.first].end();)
    {
      fptr fptr_ele1 = *it;
      if (fptr_ele1->detect == TRUE) {
        //cout <<"drop11: "  <<fptr_ele1->fault_no << " " << fptr_ele1->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[fptr_ele1->to_swlist]->name << "SA" << fptr_ele1->fault_type << endl;

        if(std::find(PNF1[v.first].begin(), PNF1[v.first].end(), fptr_ele1) == PNF1[v.first].end() ) 
        {
          PNF1[v.first].push_back(fptr_ele1);
        }
        it = SF1[v.first].erase(it);
      }
      else ++it;
    }
  }


}/* end of fault_sim_a_vector */
/* fault simulate a single test vector */
void ATPG::single_SAF_simulation23(const string &vec, int& curr) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }
  nckt = sort_wlist.size();
  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) 
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); // do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed =FALSE;
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i

  /* walk through every undetected fault
   * the undetected fault list is linked by pnext_undetect */
  //cout << "START"<< endl;
  for (int index = 0; index < PNF.size(); index ++) {
    f = PNF[index];
    //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << endl;

    f->detect = FALSE;
    if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */

    if (f->fault_type != sort_wlist[f->to_swlist]->value) {

      /* if f is a primary output or is directly connected to an primary output
       * the fault is detected */
      if ((f->node->type == OUTPUT) ||
          (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
        f->detected_time++; /// BONUS
      } else {

        /* if f is an gate output fault */
        if (f->io == GO) {

          /* if this wire is not yet marked as faulty, mark the wire as faulty
           * and insert the corresponding wire to the list of faulty wires. */
          if (!(sort_wlist[f->to_swlist]->is_faulty())) {
            sort_wlist[f->to_swlist]->set_faulty();
            wlist_faulty.push_front(sort_wlist[f->to_swlist]);
          }

          /* add the fault to the simulated fault list and inject the fault */
          simulated_fault_list[num_of_fault] = f;
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);

          /* mark the wire as having a fault injected
           * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
            pos_n->owire.front()->set_scheduled();
          }

          /* increment the number of simulated faults in this packet */
          num_of_fault++;
          /* start_wire_index keeps track of the smallest level of fault in this packet.
           * this saves simulation time.  */
          start_wire_index = min(start_wire_index, f->to_swlist);
        }  // if gate output fault

          /* the fault is a gate input fault */
        else {

          /* if the fault is propagated, set faulty_wire equal to the faulty wire.
           * faulty_wire is the gate output of f.  */
          faulty_wire = get_faulty_wire(f, fault_type);
          if (faulty_wire != nullptr) {

            /* if the faulty_wire is a primary output, it is detected */
            if (faulty_wire->is_output()) {
              f->detected_time++; /// BONUS
              //cout << sort_wlist[f->to_swlist]->name << " " <<f->node->name  <<endl;
              if(f->detected_time >= detected_num){
                f->detect = TRUE;
              }
            } else {
              /* if faulty_wire is not already marked as faulty, mark it as faulty
               * and add the wire to the list of faulty wires. */
              if (!(faulty_wire->is_faulty())) {
                faulty_wire->set_faulty();
                wlist_faulty.push_front(faulty_wire);
              }

              /* add the fault to the simulated list and inject it */
              simulated_fault_list[num_of_fault] = f;
              inject_fault_value(faulty_wire, num_of_fault, fault_type);

              /* mark the faulty_wire as having a fault injected
               *  and schedule the outputs of this gate */
              faulty_wire->set_fault_injected();
              for (auto pos_n : faulty_wire->onode) {
                pos_n->owire.front()->set_scheduled();
              }

              num_of_fault++;
              start_wire_index = min(start_wire_index, f->to_swlist);
            }
          }
        }
      } // if  gate input fault
    } // if fault is active


    /*
     * fault simulation of a packet
     */

    /* if this packet is full (16 faults)
     * or there is no more undetected faults remaining (pos points to the final element of flist_undetect),
     * do the fault simulation */
    if ((num_of_fault == num_of_faults_in_parallel) || (index == PNF.size() - 1)) {

      /* starting with start_wire_index, evaulate all scheduled wires
       * start_wire_index helps to save time. */
      for (i = start_wire_index; i < nckt; i++) {
        if (sort_wlist[i]->is_scheduled()) {
          sort_wlist[i]->remove_scheduled();
          fault_sim_evaluate(sort_wlist[i]);
        }
      } /* event evaluations end here */

      /* pop out all faulty wires from the wlist_faulty
        * if PO's value is different from good PO's value, and it is not unknown
        * then the fault is detected.
        *
        * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
      */
      //for (i = start_wire_index; i < nckt; i++) printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);

      while (!wlist_faulty.empty()) {
        w = wlist_faulty.front();
        wlist_faulty.pop_front();
        w->remove_faulty();
        w->remove_fault_injected();
        w->set_fault_free();
        /*TODO*/
        /*
         * After simulation is done,if wire is_output(), we should compare good value(wire_value_g) and faulty value(wire_value_f). 
         * If these two values are different and they are not unknown, then the fault is detected.  We should update the simulated_fault_list.  Set detect to true if they are different.
         * Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
         * After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g).
         */
        // 1.1. if wire is_output()
        if(w->is_output()){
          // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
          int detected = w->wire_value_g ^ w->wire_value_f;
          for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for each faults
            // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
            // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
            if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              //cout << "!23 "<< endl;
              if(mismatching_output.find(w->name) != mismatching_output.end()) simulated_fault_list[f_idx]->detect = TRUE;
            }
          }
        }
        // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
        w->wire_value_f = w->wire_value_g;
        /*end of TODO*/

      } // pop out all faulty wires
      num_of_fault = 0;  // reset the counter of faults in a packet
      start_wire_index = 10000;  //reset this index to a very large value.
    } // end fault sim of a packet
  } // end loop. for f = flist
  
  PNF.erase(
    std::remove_if(PNF.begin(), PNF.end(),
        [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          
          //cout <<"drop2: "  <<fptr_ele->fault_no << " " << fptr_ele->node->name << ":" << (fptr_ele->io?"O":"I")<< " "  << sort_wlist[fptr_ele->to_swlist]->name << "SA" << fptr_ele->fault_type << endl;
          curr++;
          SF.push_back(fptr_ele);
          return true;
        } else {
          return false;
        }
      }),
    PNF.end());
  /* fault dropping  */
  for(auto v : PNF1)
  {
    for (auto it = PNF1[v.first].begin(); it != PNF1[v.first].end();)
    {
      fptr fptr_ele1 = *it;
      if (fptr_ele1->detect == TRUE) {
        //cout <<"drop11: "  <<fptr_ele1->fault_no << " " << fptr_ele1->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[fptr_ele1->to_swlist]->name << "SA" << fptr_ele1->fault_type << endl;

        if(std::find(SF1[v.first].begin(), SF1[v.first].end(), fptr_ele1) == SF1[v.first].end() ) 
        {
          SF1[v.first].push_back(fptr_ele1);
        }
        it = PNF1[v.first].erase(it);
      }
      else ++it;
    }
  }
  
}/* end of fault_sim_a_vector */


/* fault simulate a single test vector */
void ATPG::single_SAF_simulation34(const string &vec, int& curr) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  nckt = sort_wlist.size();
  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed = FALSE;
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i

  /* walk through every undetected fault
   * the undetected fault list is linked by pnext_undetect */
  for (int index = 0; index < SF.size(); index ++) {
    f = SF[index];
    f->detect = FALSE;
    if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */
    
    /* consider only active (aka. excited) fault
     * (sa1 with correct output of 0 or sa0 with correct output of 1) */
    if (f->fault_type != sort_wlist[f->to_swlist]->value) {

      /* if f is a primary output or is directly connected to an primary output
       * the fault is detected */
      if ((f->node->type == OUTPUT) ||
          (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
      } else {

        /* if f is an gate output fault */
        if (f->io == GO) {

          /* if this wire is not yet marked as faulty, mark the wire as faulty
           * and insert the corresponding wire to the list of faulty wires. */
          if (!(sort_wlist[f->to_swlist]->is_faulty())) {
            sort_wlist[f->to_swlist]->set_faulty();
            wlist_faulty.push_front(sort_wlist[f->to_swlist]);
          }

          /* add the fault to the simulated fault list and inject the fault */
          simulated_fault_list[num_of_fault] = f;
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);

          /* mark the wire as having a fault injected
           * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
            pos_n->owire.front()->set_scheduled();
          }

          /* increment the number of simulated faults in this packet */
          num_of_fault++;
          /* start_wire_index keeps track of the smallest level of fault in this packet.
           * this saves simulation time.  */
          start_wire_index = min(start_wire_index, f->to_swlist);
        }  // if gate output fault

          /* the fault is a gate input fault */
        else {

          /* if the fault is propagated, set faulty_wire equal to the faulty wire.
           * faulty_wire is the gate output of f.  */
          faulty_wire = get_faulty_wire(f, fault_type);
          if (faulty_wire != nullptr) {

            /* if the faulty_wire is a primary output, it is detected */
            if (faulty_wire->is_output() ) {
                if(mismatching_output.find(faulty_wire->name) != mismatching_output.end())f->detect = TRUE;
            } else {
              /* if faulty_wire is not already marked as faulty, mark it as faulty
               * and add the wire to the list of faulty wires. */
              if (!(faulty_wire->is_faulty())) {
                faulty_wire->set_faulty();
                wlist_faulty.push_front(faulty_wire);
              }

              /* add the fault to the simulated list and inject it */
              simulated_fault_list[num_of_fault] = f;
              inject_fault_value(faulty_wire, num_of_fault, fault_type);

              /* mark the faulty_wire as having a fault injected
               *  and schedule the outputs of this gate */
              faulty_wire->set_fault_injected();
              for (auto pos_n : faulty_wire->onode) {
                pos_n->owire.front()->set_scheduled();
              }

              num_of_fault++;
              start_wire_index = min(start_wire_index, f->to_swlist);
            }
          }
        }
      } // if  gate input fault
    } // if fault is active


    /*
     * fault simulation of a packet
     */

    /* if this packet is full (16 faults)
     * or there is no more undetected faults remaining (pos points to the final element of flist_undetect),
     * do the fault simulation */
    if ((num_of_fault == num_of_faults_in_parallel) || (index == SF.size() - 1)) {

      /* starting with start_wire_index, evaulate all scheduled wires
       * start_wire_index helps to save time. */
      for (i = start_wire_index; i < nckt; i++) {
        if (sort_wlist[i]->is_scheduled()) {
          sort_wlist[i]->remove_scheduled();
          fault_sim_evaluate(sort_wlist[i]);
        }
      } /* event evaluations end here */
      //for (i = start_wire_index; i < nckt; i++) printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);

      /* pop out all faulty wires from the wlist_faulty
        * if PO's value is different from good PO's value, and it is not unknown
        * then the fault is detected.
        *
        * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
      */
      while (!wlist_faulty.empty()) {
        w = wlist_faulty.front();
        wlist_faulty.pop_front();
        w->remove_faulty();
        w->remove_fault_injected();
        w->set_fault_free();
        /*TODO*/
        /*
         * After simulation is done,if wire is_output(), we should compare good value(wire_value_g) and faulty value(wire_value_f). 
         * If these two values are different and they are not unknown, then the fault is detected.  We should update the simulated_fault_list.  Set detect to true if they are different.
         * Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
         * After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g).
         */
        // 1.1. if wire is_output()
        if(w->is_output()){
          // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
          int detected = w->wire_value_g ^ w->wire_value_f;
          for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for each faults
            // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
            // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
            if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              
              if(mismatching_output.find(w->name) != mismatching_output.end()){
                  simulated_fault_list[f_idx]->detect = TRUE;
              } 
            }
          }
        }
        // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
        w->wire_value_f = w->wire_value_g;
        /*end of TODO*/

      } // pop out all faulty wires
      num_of_fault = 0;  // reset the counter of faults in a packet
      start_wire_index = 10000;  //reset this index to a very large value.
    } // end fault sim of a packet
  } // end loop. for f = flist
  SF.erase(
    std::remove_if(SF.begin(), SF.end(),
        [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          //cout <<"drop6: "  <<fptr_ele->fault_no << " " << fptr_ele->node->name << ":" << (fptr_ele->io?"O":"I")<< " "  << sort_wlist[fptr_ele->to_swlist]->name << "SA" << fptr_ele->fault_type << endl;
          curr++;
          PNF.push_back(fptr_ele);
          //cout << sort_wlist[fptr_ele->to_swlist]->name << " " <<fptr_ele->detect <<  endl;

          return true;
        } else if (fptr_ele->detect == FALSE){
          return false;
        }
      }),
    SF.end());
  
  for(auto v : SF1)
  {
    for (auto it = SF1[v.first].begin(); it != SF1[v.first].end();)
    {
      fptr fptr_ele1 = *it;
      if (fptr_ele1->detect == TRUE) {
        //cout <<"drop11: "  <<fptr_ele1->fault_no << " " << fptr_ele1->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[fptr_ele1->to_swlist]->name << "SA" << fptr_ele1->fault_type << endl;

        if(std::find(PNF1[v.first].begin(), PNF1[v.first].end(), fptr_ele1) == PNF1[v.first].end() ) 
        {
          PNF1[v.first].push_back(fptr_ele1);
        }
        it = SF1[v.first].erase(it);
      }
      else ++it;
    }
  }
  /* fault dropping  */

}/* end of fault_sim_a_vector */


/* multiple SAF simulate a single test vector */
bool ATPG::multiple_SAF_simulation2(const string &vec, fptr examined_fault) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;
  bool flag;
  num_of_fault = 0; // counts the number of faults in a packet


  /* Keep track of the minimum wire index of 16 faults in a packet.
   * the start_wire_index is used to keep track of the
   * first gate that needs to be evaluated.
   * This reduces unnecessary check of scheduled events.*/
  start_wire_index = 10000;

  /* for every input, set its value to the current vector value */
  for (i = 0; i < cktin.size(); i++) {
    cktin[i]->value = ctoi(vec[i]);
  }

  /* initialize the circuit - mark all inputs as changed and all other
   * nodes as unknown (2) */
  nckt = sort_wlist.size();
  for (i = 0; i < nckt; i++) {
    if (i < cktin.size()) {
      sort_wlist[i]->set_changed();
    } else {
      sort_wlist[i]->value = U;
    }
  }

  sim(); /* do a fault-free simulation, see sim.cpp */
  if (debug) { display_io(); }

  /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
   * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
  for (i = 0; i < nckt; i++) {
    sort_wlist[i]->remove_faulty();
    sort_wlist[i]->remove_fault_injected();
    sort_wlist[i]->fixed = FALSE;
    
    switch (sort_wlist[i]->value) {
      case 1:
        sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
        sort_wlist[i]->wire_value_f = ALL_ONE;
        break;
      case 2:
        sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
        sort_wlist[i]->wire_value_f = 0x55555555;
        break;
      case 0:
        sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
        sort_wlist[i]->wire_value_f = ALL_ZERO;
        break;
    }
  } // for i
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    f->detect = false;
    if(examined_fault == f) continue;
    if ((f->node->type == OUTPUT) ||
        (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
      if (!(sort_wlist[f->to_swlist]->is_faulty())) {
        sort_wlist[f->to_swlist]->set_faulty();
        wlist_faulty.push_front(sort_wlist[f->to_swlist]);
      }
      if (sort_wlist[f->to_swlist]->fixed == FALSE){ 
        inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
        /* mark the wire as having a fault injected
          * and schedule the outputs of this gate */
        sort_wlist[f->to_swlist]->set_fault_injected();
        sort_wlist[f->to_swlist]->fixed = TRUE;
        start_wire_index = min(start_wire_index, f->to_swlist);
      } 
      else
        {
          //f->detect =TRUE;
        }
    }
    else 
    {
      if (f->io == GO) {
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f << " " << f->fault_type  << " fixed? " << sort_wlist[f->to_swlist]->fixed<<endl;
        if (!(sort_wlist[f->to_swlist]->is_faulty())) {
          sort_wlist[f->to_swlist]->set_faulty();
          wlist_faulty.push_front(sort_wlist[f->to_swlist]);
        }
        if (sort_wlist[f->to_swlist]->fixed == FALSE){
          //cout << " inject "<< endl; 
          inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
          /* mark the wire as having a fault injected
            * and schedule the outputs of this gate */
          sort_wlist[f->to_swlist]->set_fault_injected();
          sort_wlist[f->to_swlist]->fixed = TRUE;
          for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
        else
        {
          //f->detect =TRUE;
        }
        //cout << sort_wlist[f->to_swlist]->name <<  sort_wlist[f->to_swlist]->wire_value_f <<endl;
      }
      else
      { 
        continue;
      }
    }
  }
  for (i = 0; i < 0; i++) {
    printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
  }
  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  } /* event evaluations end here */

  
  start_wire_index = 10000;
  for (int index = 0; index < SF.size(); index ++) {
    f =  SF[index];
    if(examined_fault == f) continue;
    if ((f->node->type == OUTPUT) ||
        (f->io == GO )) {
        continue;
    }
    else // Gate input
    {
      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f1 = wire->wire_value_f;
            inject_fault_value(wire, num_of_fault, ff->fault_type);
            
          }
        }
      }
      //cout << sort_wlist[ f->to_swlist]->name.c_str() << sort_wlist[ f->to_swlist]->wire_value_f1 << " " << sort_wlist[ f->to_swlist]->wire_value_f << endl;

      sort_wlist[ f->to_swlist]->wire_value_f1 = sort_wlist[ f->to_swlist]->wire_value_f;

      faulty_wire = get_faulty_wire(f, fault_type);
      //faulty_wire = get_faulty_wire(f, fault_type);
      if (faulty_wire != nullptr){
        sort_wlist[ f->to_swlist]->wire_value_f1 = sort_wlist[ f->to_swlist]->wire_value_f;
        inject_fault_value( sort_wlist[ f->to_swlist], num_of_fault, f->fault_type);
        faulty_wire = f->node->owire.front();
        
        if(faulty_wire->fixed == FALSE) fault_sim_evaluate(faulty_wire);
        else {
          //cout <<  sort_wlist[ f->to_swlist]->name << " " << faulty_wire->name << faulty_wire->fixed << endl;
          //f->detect = TRUE;
        }
        sort_wlist[ f->to_swlist]->wire_value_f = sort_wlist[ f->to_swlist]->wire_value_f1;
      }
      else
      {
        //cout <<  sort_wlist[ f->to_swlist]->name << " null " << endl;
      }
      
      sort_wlist[ f->to_swlist]->wire_value_f = sort_wlist[ f->to_swlist]->wire_value_f1;
      //cout << sort_wlist[ f->to_swlist]->name.c_str() << sort_wlist[ f->to_swlist]->wire_value_f1 << " " << sort_wlist[ f->to_swlist]->wire_value_f << endl;

      for (auto wire : f->node->iwire)
      {
        //printf("%s %x %s %x ",   wire->name .c_str(), wire->wire_value_f,   sort_wlist[ f->to_swlist]->name .c_str(), sort_wlist[ f->to_swlist]->wire_value_f);
        for(auto ff : SF1[wire->name])
        {
          if ((ff->io == 0) && (ff->node->name == f->node->name) && (wire->name != sort_wlist[ f->to_swlist]->name)){
            wire->wire_value_f = wire->wire_value_f1;
          }
        }
      }
      if (faulty_wire != nullptr ) {
        /* if the faulty_wire is a primary output, it is detected */
        if (faulty_wire->fixed == FALSE ) {
          if (!(faulty_wire->is_faulty())) {
            faulty_wire->set_faulty();
            wlist_faulty.push_front(faulty_wire);
          }
          /* add the fault to the simulated list and inject it */
          //inject_fault_value(faulty_wire, num_of_fault, fault_type);
          /* mark the faulty_wire as having a fault injected
          *  and schedule the outputs of this gate */
          faulty_wire->set_fault_injected();
          faulty_wire->fixed == TRUE;
          if(!faulty_wire->is_output()){
            for (auto pos_n : faulty_wire->onode) {
              if(pos_n->owire.front()->fixed == FALSE)pos_n->owire.front()->set_scheduled();
            }
          }
          start_wire_index = min(start_wire_index, f->to_swlist);
        }
      }
    }
  }

  for (i = start_wire_index; i < nckt; i++) {
    if (sort_wlist[i]->is_scheduled()) {
      sort_wlist[i]->remove_scheduled();
      fault_sim_evaluate(sort_wlist[i]);
    }
  }

  /* pop out all faulty wires from the wlist_faulty
    * if PO's value is different from good PO's value, and it is not unknown
    * then the fault is detected.
    *
    * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
  */
  for (i = 0; i < 0; i++) {
    printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
  }
  for (auto ww : cktout) {
    //printf("%s %x %x\n",  sort_wlist[i]->name.c_str(), sort_wlist[i]->wire_value_f, sort_wlist[i]->wire_value_g);
    //cout << sort_wlist[i]->name << " " << sort_wlist[i]->wire_value_f&0x03 << endl;
    w = ww;
    if(w->is_output()){
      // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
      int detected = w->wire_value_g ^ w->wire_value_f;
      //cout << w->name << vec << " "<< detected << endl;
      for(int f_idx=0; f_idx<1; ++f_idx){ // for each faults
        // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
        // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
        if(((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
          // If sim result is different from observed value -> record mismatching output
          //cout <<w->name << w->wire_value_f << w->wire_value_g <<endl;
          if(((Mask[f_idx] & detected) != 0)&&(tr_unexamined1[vec].find(w->name) != tr_unexamined1[vec].end())) {
            //cout << w->name << vec << " "<< detected << endl;
            flag = true;
          }
          else if (((Mask[f_idx] & detected) == 0)&&(tr_unexamined1[vec].find(w->name) == tr_unexamined1[vec].end())){
            //cout << vec<< w->name  << endl;
            flag =  true;
          }
          else
          {
            
            return false;
          }
        }

      }
    }
  }

  return flag;
  while (!wlist_faulty.empty()) {
    w = wlist_faulty.front();
    wlist_faulty.pop_front();
    w->remove_faulty();
    w->remove_fault_injected();
    w->set_fault_free();
    // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
    w->wire_value_f = w->wire_value_g;


  } // pop out all faulty wires

  //for (auto m :mismatching_output ) cout << m << endl;
  num_of_fault = 0;  // reset the counter of faults in a packet
  start_wire_index = 10000;  //reset this index to a very large value.
}/* end of fault_sim_a_vector */

void ATPG::diag1(){
  string  wire_name;
  int observed, expected;
  trptr fc;
  trptr temp;
  vector<string> test;
  test.push_back("01000");
   //test.push_back("11011");
    //test.push_back("01100");
  //reverse(SF.begin(), SF.end());
  for (auto vec : vectors) {
    //vec = tc.first;
    
    wptr w, faulty_wire;
    /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
    fptr simulated_fault_list[num_of_faults_in_parallel];
    fptr f;
    int fault_type;
    int i, start_wire_index, nckt;
    int num_of_fault;
    num_of_fault = 0; // counts the number of faults in a packet

    /* Keep track of the minimum wire index of 16 faults in a packet.
    * the start_wire_index is used to keep track of the
    * first gate that needs to be evaluated.
    * This reduces unnecessary check of scheduled events.*/
    start_wire_index = 10000;
  
    /* for every input, set its value to the current vector value */
    for (i = 0; i < cktin.size(); i++) {
      cktin[i]->value = ctoi(vec[i]);
    }

    /* initialize the circuit - mark all inputs as changed and all other
    * nodes as unknown (2) */
    nckt = sort_wlist.size();
    for (i = 0; i < nckt; i++) {
      if (i < cktin.size()) {
        sort_wlist[i]->set_changed();
      } else {
        sort_wlist[i]->value = U;
      }
    }

    sim(); /* do a fault-free simulation, see sim.cpp */
    /* expand the fault-free value into 32 bits (00 = logic zero, 11 = logic one, 01 = unknown)
    * and store it in wire_value_g (good value) and wire_value_f (faulty value)*/
    for (i = 0; i < nckt; i++) {
      sort_wlist[i]->remove_scheduled();
      sort_wlist[i]->remove_faulty();
      sort_wlist[i]->remove_fault_injected();
      switch (sort_wlist[i]->value) {
        
        case 1:
          sort_wlist[i]->wire_value_g = ALL_ONE;  // 11 represents logic one
          sort_wlist[i]->wire_value_f = ALL_ONE;
          break;
        case 2:
          sort_wlist[i]->wire_value_g = 0x55555555; // 01 represents unknown
          sort_wlist[i]->wire_value_f = 0x55555555;
          break;
        case 0:
          sort_wlist[i]->wire_value_g = ALL_ZERO; // 00 represents logic zero
          sort_wlist[i]->wire_value_f = ALL_ZERO;
          break;
      }
    }
    for (int index = 0; index < SF.size(); index ++){
      f = SF[index];
      //cout << vec << " fault no: " << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << " tfsf: " << f->tfsf << " tpsf: " << f->tpsf << endl;
      
      if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */
      /* consider only active (aka. excited) fault
       * (sa1 with correct output of 0 or sa0 with correct output of 1) */
      // if (...) // TODO: Step 1. Structural Backtracing: Check whether pos is in all failing fanin cone
      // if (...) // TODO: Step 2. Parity Check: Check whether vp=f
      if ( f->fault_type != sort_wlist[f->to_swlist]->value) { // Step 3. Excitation Condition Check: Check the sa-v != n
        /* if f is a primary output or is directly connected to an primary output
         * the fault is detected */
        if ((f->node->type == OUTPUT) ||
            (f->io == GO && sort_wlist[f->to_swlist]->is_output())) {
          
          if(!tr_unexamined1[vec].empty())
          {
            if((tr_unexamined1[vec].find(sort_wlist[f->to_swlist]->name) != tr_unexamined1[vec].end())) { 
              f->tfsf++;
              f->tfsp--;
            }
            else {f->tpsf++;}
            // no fault is supposed to occur
          } else {f->tpsf++;}
        
        } else {
          /* if f is an gate output fault */
          if (f->io == GO) {
            /* if this wire is not yet marked as faulty, mark the wire as faulty
             * and insert the corresponding wire to the list of faulty wires. */
            if (!(sort_wlist[f->to_swlist]->is_faulty())) {
              sort_wlist[f->to_swlist]->set_faulty();
              wlist_faulty.push_front(sort_wlist[f->to_swlist]);
            }
  
            /* add the fault to the simulated fault list and inject the fault */
            simulated_fault_list[num_of_fault] = f;
            inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
  
            /* mark the wire as having a fault injected
             * and schedule the outputs of this gate */
            sort_wlist[f->to_swlist]->set_fault_injected();
            for (auto pos_n : sort_wlist[f->to_swlist]->onode) {
              pos_n->owire.front()->set_scheduled();
            }
  
            /* increment the number of simulated faults in this packet */
            num_of_fault++;
            /* start_wire_index keeps track of the smallest level of fault in this packet.
             * this saves simulation time.  */
            start_wire_index = min(start_wire_index, f->to_swlist);
          }  // if gate output fault
  
          /* the fault is a gate input fault */
          else {
            /* if the fault is propagated, set faulty_wire equal to the faulty wire.
             * faulty_wire is the gate output of f.  */
            //cout << "123 " << endl;
            faulty_wire = get_faulty_wire(f, fault_type);
            if (faulty_wire != nullptr) {
              /* if the faulty_wire is a primary output, it is detected */
              if (faulty_wire->is_output()) {
                if(!tr_unexamined1[vec].empty())
                {
                  if((tr_unexamined1[vec].find(faulty_wire->name) != tr_unexamined1[vec].end())) { 
                    f->tfsf++;
                    f->tfsp--;
                  }
                  else {f->tpsf++;}
                  // no fault is supposed to occur
                } else {f->tpsf++;}
              } else {
                /* if faulty_wire is not already marked as faulty, mark it as faulty
                 * and add the wire to the list of faulty wires. */
                if (!(faulty_wire->is_faulty())) {
                  faulty_wire->set_faulty();
                  wlist_faulty.push_front(faulty_wire);
                }
  
                /* add the fault to the simulated list and inject it */
                simulated_fault_list[num_of_fault] = f;
                inject_fault_value(faulty_wire, num_of_fault, fault_type);
  
                /* mark the faulty_wire as having a fault injected
                 *  and schedule the outputs of this gate */
                faulty_wire->set_fault_injected();
                for (auto pos_n : faulty_wire->onode) {
                  pos_n->owire.front()->set_scheduled();
                }
  
                num_of_fault++;
                start_wire_index = min(start_wire_index, f->to_swlist);
              }
            }
          }
        } // if  gate input fault
      } // if fault is active
      else f->tpsp++;
      /*
       * fault simulation of a packet
       */
  
      /* if this packet is full (16 faults)
       * or there is no more undetected faults remaining (pos points to the final element of flist_undetect),
       * do the fault simulation */
      if ((num_of_fault == num_of_faults_in_parallel) || (index == SF.size()-1)) {
  
        /* starting with start_wire_index, evaulate all scheduled wires
         * start_wire_index helps to save time. */
        //cout << nckt << endl;
        for (i = start_wire_index; i < nckt; i++) {
          
          if (sort_wlist[i]->is_scheduled()) {
            //cout << sort_wlist[i]->name << sort_wlist[i]->inode.front()->name << endl;
            sort_wlist[i]->remove_scheduled();
            fault_sim_evaluate(sort_wlist[i]);
          }
        } /* event evaluations end here */
  
        /* pop out all faulty wires from the wlist_faulty
          * if PO's value is different from good PO's value, and it is not unknown
          * then the fault is detected.
          *
          * IMPORTANT! remember to reset the wires' faulty values back to fault-free values.
        */
        while (!wlist_faulty.empty()) {
          w = wlist_faulty.front();
          wlist_faulty.pop_front();
          w->remove_faulty();
          w->remove_fault_injected();
          w->set_fault_free();
          /*TODO*/
          /*
           * After simulation is done,if wire is_output(), we should compare good value(wire_value_g) and faulty value(wire_value_f). 
           * If these two values are different and they are not unknown, then the fault is detected.  We should update the simulated_fault_list.  Set detect to true if they are different.
           * Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
           * After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g).
           */
          // 1.1. if wire is_output()
          if(w->is_output()){
            //cout << w->name << endl;
            
              
            // 1.2. we should compare good value(wire_value_g) and faulty value(wire_value_f)
            int detected = w->wire_value_g ^ w->wire_value_f;
            for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for every undetected fault
  
              // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
              // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
                if(!tr_unexamined1[vec].empty())
                {
                  if((tr_unexamined1[vec].find(w->name) != tr_unexamined1[vec].end())) { 
                    simulated_fault_list[f_idx]->tfsf++;
                    simulated_fault_list[f_idx]->tfsp--;
                  }
                  else {simulated_fault_list[f_idx]->tpsf++;}
                  // no fault is supposed to occur
                } else {simulated_fault_list[f_idx]->tpsf++;}
              }
              //cout << endl;
            }
          }
          // 4. After that, don't forget to reset faulty values (wire_value_f) to their fault-free values (wire_value_g)
          w->wire_value_f = w->wire_value_g;
          /*end of TODO*/
        } // pop out all faulty wires
        num_of_fault = 0;  // reset the counter of faults in a packet
        start_wire_index = 10000;  //reset this index to a very large value.
      } // end fault sim of a packet
    }
  
  }

  for(auto ff : SF)
  {
    ff->score = ((double) ff->tfsf*10)/((double) (ff->tfsf*10 + ff->tfsp*10 + ff->tpsf))*100;
    
  }

}