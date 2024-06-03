/**********************************************************************/
/*           Parallel-Fault Event-Driven Fault Simulator              */
/*                                                                    */
/*           Author: Bing-Chen (Benson) Wu                            */
/*           last update : 02/19/2024                                 */
/**********************************************************************/

#include "atpg.h"

/* pack 16 faults into one packet.  simulate 16 faults together. */
#define num_of_faults_in_parallel  16

/* The faulty_wire contains a list of wires that 
 * change their values in the fault simulation for a particular packet.
 * (that is wire_value_g != wire_value_f) 
 * Note that the wire themselves are not necessarily a fault site.
 * The list is linked by the pnext pointers */

/* fault simulate a set of test vectors */
void ATPG::fault_simulate_vectors(int &total_detect_num) {
  int i;
  int current_detect_num = 0;

  /* for every vector */
  for (i = vectors.size() - 1; i >= 0; i--) {
    // cout << "fault_sim a vector(" << vectors[i] << "):\n";
    fault_sim_a_vector(vectors[i], current_detect_num);
    total_detect_num += current_detect_num;
    fprintf(stdout, "vector[%d] detects %d faults (%d)\n", i, current_detect_num, total_detect_num);
  }
}// fault_simulate_vectors

/* fault simulate a single test vector */
void ATPG::fault_sim_a_vector(const string &vec, int &num_of_current_detect) {
  wptr w, faulty_wire;
  /* array of 16 fptrs, which points to the 16 faults in a simulation packet  */
  fptr simulated_fault_list[num_of_faults_in_parallel];
  fptr f;
  int fault_type;
  int i, start_wire_index, nckt;
  int num_of_fault;

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
        f->detected_time++; /// BONUS
        if(f->detected_time >= detected_num){
          f->detect = TRUE;
        }
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
          for(int f_idx=0; f_idx<num_of_fault; ++f_idx){ // for each faults
            // 2. If these two values are different and they are not unknown, then the fault is detected. Since we use two-bit logic to simulate circuit, you can use Mask[] to perform bit-wise operation to get value of a specific bit.
            // if((Mask[f_idx] & detected != 0) && (w->wire_value_g & Mask[f_idx] != Unknown[f_idx]) && (w->wire_value_f & Mask[f_idx] != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
            if(((Mask[f_idx] & detected) != 0) && ((w->wire_value_g & Mask[f_idx]) != Unknown[f_idx]) && ((w->wire_value_f & Mask[f_idx]) != Unknown[f_idx])){ // the bits are faults and are detected. And wire_value_g and wire_value_g are not unknown.
              simulated_fault_list[f_idx]->detect = TRUE;
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

  /* fault dropping  */
  flist_undetect.remove_if(
      [&](const fptr fptr_ele) {
        if (fptr_ele->detect == TRUE) {
          num_of_current_detect += fptr_ele->eqv_fault_num;
          return true;
        } else {
          return false;
        }
      });
}/* end of fault_sim_a_vector */

/* evaluate wire w 
 * 1. update w->wire_value_f 
 * 2. schedule new events if value2 != value1 */
void ATPG::fault_sim_evaluate(const wptr w) {
  unsigned int new_value;
  nptr n;
  int i, nin, nout;

  n = w->inode.front();
  nin = n->iwire.size();
  switch (n->type) {
    /*break a multiple-input gate into multiple two-input gates */
    case AND:
    case BUF:
    case NAND:
      new_value = ALL_ONE;
      for (i = 0; i < nin; i++) {
        new_value &= n->iwire[i]->wire_value_f;
      }
      if (n->type == NAND) {
        new_value = PINV(new_value);  // PINV is for three-valued inversion
      }
      break;

    case OR:
    case NOR:
      new_value = ALL_ZERO;
      for (i = 0; i < nin; i++) {
        new_value |= n->iwire[i]->wire_value_f;
      }
      if (n->type == NOR) {
        new_value = PINV(new_value);
      }
      break;

    case NOT:
      new_value = PINV(n->iwire.front()->wire_value_f);
      break;

    case XOR:
      new_value = PEXOR(n->iwire[0]->wire_value_f, n->iwire[1]->wire_value_f);
      break;

    case EQV:
      new_value = PEQUIV(n->iwire[0]->wire_value_f, n->iwire[1]->wire_value_f);
      break;
  }

  /* if the new_value is different than the wire_value_g (the good value),
   * save it */
  if (w->wire_value_g != new_value) {

    /* if this wire is faulty, make sure the fault remains injected */
    if (w->is_fault_injected()) {
      combine(w, new_value);
    }

    /* update wire_value_f */
    w->wire_value_f = new_value;

    /* insert wire w into the faulty_wire list */
    if (!(w->is_faulty())) {
      w->set_faulty();
      wlist_faulty.push_front(w);
    }

    /* schedule new events */
    for (i = 0, nout = w->onode.size(); i < nout; i++) {
      if (w->onode[i]->type != OUTPUT) {
        w->onode[i]->owire.front()->set_scheduled();
      }
    }
  } // if new_value is different
  // if new_value is the same as the good value, do not schedule any new event
}/* end of fault_sim_evaluate */

/* Given a gate-input fault f, check if f is propagated to the gate output.
   If so, returns the gate output wire as the faulty_wire.
   Also returns the gate output fault type. 
*/
ATPG::wptr ATPG::get_faulty_wire(const fptr f, int &fault_type) {
  int i, nin;
  bool is_faulty;

  is_faulty = true;
  nin = f->node->iwire.size();
  switch (f->node->type) {

    /* this case should not occur,
     * because we do not create fault in the NOT BUF gate input */
    case NOT:
    case BUF:
      fprintf(stdout, "something is fishy(get_faulty_net)...\n");
      break;

    case AND:
      /*TODO*/
      /* To check if this fault can propagate successfully, we need to check every gate input of this
       * AND gate. You sould avoid checking the faultly input. If any input is zero or unknown, fault 
       * f can't propagate, and you sould set (is_faulty) to false
       */
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) { // fanin that not the faulty wire
          if (f->node->iwire[i]->value == 0 || f->node->iwire[i]->value == 2) { // If any input is zero or unknown, fault f can't propagate
            is_faulty = false;
          }
        }
      }
      /*end of TODO*/
      /* AND gate input stuck-at one/zero fault is propagated to
         AND gate output stuck-at one/zero fault */
      if (f->fault_type == 0)
        fault_type = STR;
      else
        fault_type = STF;
      break;

    case NAND:
      /*TODO*/
      /* To check if this fault can propagate successfully, we need to check every gate input of this
       * NAND gate. You sould avoid checking the faultly input. You sould figure out in which condition 
       * f can't propagate, and you sould set (is_faulty) to false
       */
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) { // fanin that not the faulty wire
          if (f->node->iwire[i]->value == 0 || f->node->iwire[i]->value == 2) { // If any input is zero or unknown, fault f can't propagate
            is_faulty = false;
          }
        }
      }
      /*end of TODO*/
      if (f->fault_type == 0)
        fault_type = STF;
      else
        fault_type = STR;
      break;
    case OR:
      /*TODO*/
      /* To check if this fault can propagate successfully, we need to check every gate input of this
       * OR gate. You sould avoid checking the faultly input. You sould figure out in which condition 
       * f can't propagate, and you sould set (is_faulty) to false
       */
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) { // fanin that not the faulty wire
          if (f->node->iwire[i]->value == 1 || f->node->iwire[i]->value == 2) { // If any input is one or unknown, fault f can't propagate
            is_faulty = false;
          }
        }
      }
      /*end of TODO*/
      if (f->fault_type == 0)
        fault_type = STR;
      else
        fault_type = STF;
      break;
    case NOR:
      /*TODO*/
      /* To check if this fault can propagate successfully, we need to check every gate input of this
       * NOR gate. You sould avoid checking the faultly input. You sould figure out in which condition 
       * f can't propagate, and you sould set (is_faulty) to false
       */
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) { // fanin that not the faulty wire
          if (f->node->iwire[i]->value == 1 || f->node->iwire[i]->value == 2) { // If any input is one or unknown, fault f can't propagate
            is_faulty = false;
          }
        }
      }
      /*end of TODO*/
      if (f->fault_type == 0)
        fault_type = STF;
      else
        fault_type = STR;
      break;
    case XOR:
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) {
          if (f->node->iwire[i]->value == 0) {
            fault_type = f->fault_type;
          } else {
            fault_type = f->fault_type ^ 1;
          }
        }
      }
      break;
    case EQV:
      for (i = 0; i < nin; i++) {
        if (f->node->iwire[i] != sort_wlist[f->to_swlist]) {
          if (f->node->iwire[i]->value == 0) {
            fault_type = f->fault_type ^ 1;
          } else {
            fault_type = f->fault_type;
          }
        }
      }
      break;
  }
  if (is_faulty) {
    return (f->node->owire.front());
  }
  return (nullptr);
}/* end of get_faulty_wire */

/* This function injects a fault 
 * 32 bits(16 faults) in a word, 
 * Mask[bit position]='11' is the bit position of the fault
 * for example, Mask[0]=0x00000003 means that the right most two bits are the position of the fault.  
 */
void ATPG::inject_fault_value(const wptr faulty_wire, const int &bit_position, const int &fault_type) {
  /*TODO*/
  //Hint: use mask to inject fault to the right position
  /* Use mask[] to perform bit operation to inject fault (STUCK1 or STUCK0) to the right position
   * Call inject_fault_at() to set the fault_flag of the injected bit position.
   */
  // 1. Use mask[] to perform bit operation to inject fault (STUCK1 or STUCK0) to the right position
  if(fault_type == STUCK0){ // SA0: use AND GATE with input 0 to simulate as stuck at '0'. (See ch5 p.21)
    faulty_wire->wire_value_f = faulty_wire->wire_value_f & ~Mask[bit_position]; 
  } else if(fault_type == STUCK1){ // SA0: use OR GATE with input 1 to simulate as stuck at '1'. (See ch5 p.21)
    faulty_wire->wire_value_f = faulty_wire->wire_value_f | Mask[bit_position];
  } else{
    cout << "[Error] - Niether SA0 nor SA1.\n";
  }
  // 2. Call inject_fault_at() to set the fault_flag of the injected bit position.
  faulty_wire->inject_fault_at(bit_position);
  /*end of TODO*/
}/* end of inject_fault_value */

/* For each fault in this packet, check
 * if wire w itself is the fault site, do not change its wire_value_f.
 * (because the wire_value_f was already decided by the inject_fault_value function) */
void ATPG::combine(const wptr w, unsigned int &new_value) {
  int i;

  for (i = 0; i < num_of_faults_in_parallel; i++) {
    if (w->has_fault_at(i)) {
      new_value &= ~Mask[i];
      new_value |= (w->wire_value_f & Mask[i]);
    }
  }
} /* end of combine */

/* for three-valued logic inversion
 * Swap the odd bits and the even bits,
 * then do a inversion (using XOR).
 * The purpose of swapping is to keep the unknown u=01 unchanged after inversion */
unsigned int ATPG::PINV(const unsigned int &value) {
  return ((((value & 0x55555555) << 1) ^ 0xaaaaaaaa) |
          (((value & 0xaaaaaaaa) >> 1) ^ 0x55555555));
}/* end of PINV */


unsigned int ATPG::PEXOR(const unsigned int &value1, const unsigned int &value2) {
  return ((value1 & PINV(value2)) | (PINV(value1) & value2));
}/* end of PEXOR */


unsigned int ATPG::PEQUIV(const unsigned int &value1, const unsigned int &value2) {
  return ((value1 | PINV(value2)) & (PINV(value1) | value2));
}/* end of PEQUIV */
