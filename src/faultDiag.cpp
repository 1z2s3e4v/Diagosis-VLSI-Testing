/**********************************************************************/
/*  Parallel-Fault Event-Driven Transition Delay Fault Simulator      */
/*                                                                    */
/*           Author: Tsai-Chieh Chen                                  */
/*           last update : 10/22/2018                                 */
/**********************************************************************/

#include "atpg.h"
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

/* fault simulate a set of test vectors */
void ATPG::fd_fault_simulation(int &total_detect_num) {
  int i;
  int current_detect_num = 0;

  /* for every vector */
  for (i = vectors.size() - 1; i >= 0; i--) {
    // cout << "fault_sim a vector(" << vectors[i] << "):\n";
    fd_fault_sim_a_vector(vectors[i], current_detect_num, i); // print fail log
    total_detect_num += current_detect_num;
  }
}// fd_fault_simulation

/* fault simulate a single test vector */
void ATPG::fd_fault_sim_a_vector(const string &vec, int &num_of_current_detect, int &vec_index) {
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
  for (auto pos = flist_undetect.cbegin(); pos != flist_undetect.cend(); ++pos) {
    f = *pos;
    if (f->detect == REDUNDANT) { continue; } /* ignore redundant faults */
    /* consider only active (aka. excited) fault
     * (sa1 with correct output of 0 or sa0 with correct output of 1) */
    if (f->fault_type != sort_wlist[f->to_swlist]->value) {
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
        /* add the fault to the simulated fault list and inject the fault */
        simulated_fault_list[num_of_fault] = f;
        inject_fault_value(sort_wlist[f->to_swlist], num_of_fault, f->fault_type);
        /* mark the wire as having a fault injected
         * and schedule the outputs of this gate */
        sort_wlist[f->to_swlist]->set_fault_injected();
        /* increment the number of simulated faults in this packet */
        num_of_fault++;
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
    if ((num_of_fault == num_of_faults_in_parallel) || (next(pos, 1) == flist_undetect.cend())) {

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
          for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for every undetected fault
            // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
            // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
            if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              failed_ws.insert(w);
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

  /* fault dropping  */
  /*flist_undetect.remove_if(
      [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          num_of_current_detect += fptr_ele->eqv_fault_num;
          return true;
        } else {
          return false;
        }
      });*/
  for(wptr failed_w : failed_ws){
    // printf("vector[%d] %s %s  # T'%s'\n", vec_index, failed_w->name.c_str(), failLog_str[failed_w->wire_value_g].c_str(), vec.c_str());
    if(failed_w->wire_value_g == 1){
      printf("vector[%d] %s %s  # T'%s' \n", vec_index, failed_w->name.c_str(), failLog_str[1].c_str(), vec.c_str());
    } else if (failed_w->wire_value_g == 0){
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
    flist_undetect.push_front(f.get());     // initial undetected fault list contains all faults
    flist.push_front(move(f));             // push into the fault list
  }
  
  flist.reverse();
  flist_undetect.reverse();
  /*walk through all faults, assign fault_no one by one  */
  fault_num = 0;
  for (fptr f : flist_undetect) {
    f->fault_no = fault_num;
    fault_num++;
    // cout << f->fault_no << " "  << f->node->name << ":" << (f->io?"O":"I") <<" "  << (f->io?9:(f->index)) <<" "  << "SA" << f->fault_type << endl;
  }
  // fprintf(stdout, "#number of equivalent faults = %d\n", fault_num);
}

void ATPG::read_faillog(const string &faillog) {  
  /*
  string t, vec;
  size_t i;
  ifstream file(vetFile, std::ifstream::in); // open the input vectors' file
  if (!file) { // if the ifstream obj does not exist, fail to open the file
    fprintf(stderr, "File %s could not be opened\n", vetFile.c_str());
    exit(EXIT_FAILURE);
  }
    while (!file.eof() && !file.bad()) {
    getline(file, t); // get a line from the file
    vec.clear();
    for (char c: t) {
      if (c == 'T') continue;  // ignore "T"
      if (c == '\'') continue; // ignore "'"
      if (c == ' ') continue;  // ignore " "
      vec.push_back(c);
    //cout << "Before erase: " << t << endl;
    //cout << "After erase: " << vec << endl;
    vectors.push_back(vec); // append the vectors
    }
  }
  file.close(); // close the file
  */
}
