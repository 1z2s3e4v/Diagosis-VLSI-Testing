/**********************************************************************/
/*           Constructors of the classes defined in atpg.h            */
/*           ATPG top-level functions                                 */
/*           Author: Bing-Chen (Benson) Wu                            */
/*           last update : 02/20/2023                                 */
/**********************************************************************/

#include "atpg.h"



void ATPG::ranking()
{
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
        if (swapped == false)
            break;
    }
}

    
void ATPG::test() {
    string vec;
    int current_detect_num = 0;
    int total_detect_num = 0;
    int total_no_of_backtracks = 0;  // accumulative number of backtracks
    int current_backtracks = 0;
    int no_of_aborted_faults = 0;
    int no_of_redundant_faults = 0;
    int no_of_calls = 0;

    fptr fault_under_test;
    if (!failLog_only) fault_under_test = flist_undetect.front();
    else fault_under_test = examined_faults.front();
    
    /* stuck-at fault sim mode */
    if (fsim_only) {
        fault_simulate_vectors(total_detect_num);
        in_vector_no += vectors.size();
        display_undetect();
        fprintf(stdout, "\n");
        return;
    }// if fsim only

    /* transition fault sim mode */
    if (tdfsim_only) {
        // transition_delay_fault_simulation(total_detect_num);
        in_vector_no += vectors.size();
        display_undetect();

        printf("\n# Result:\n");
        printf("-----------------------\n");
        printf("# total transition delay faults: %d\n", num_of_tdf_fault);
        printf("# total detected faults: %d\n", total_detect_num);
        printf("# fault coverage: %lf %\n", (double) total_detect_num / (double) num_of_tdf_fault * 100);
        return;
    }// if fsim only

    /* transition fault sim mode */

    if (failLog_only) {
        fd_fault_simulation(total_detect_num);
        in_vector_no += vectors.size();
        //display_fd_undetect();
        return;
    }// if failLog only

    if (diag_only) {

        
        //display_undetect();
        diag();
        ranking();
        int count = 1;
        int group = 0;
        int max = 9999;
        for (fptr f: ranks) {
            
            if(max > f->score) {max= f->score; group++;}
            printf("No.%d %s %s %s SA%d, GroudID: %d, TFSF: %d, TPSF: %d, Score: %d     \n", count++, sort_wlist[f->to_swlist]->name.c_str(), f->node->name.c_str(), \
                                             (f->io?"GO":"GI"), f->fault_type, group, f->tfsf, f->tpsf, f->score);
            //cout << f->fault_no << " " << f->node->name << ":" << (f->io?"O":"I")<< " "  << sort_wlist[f->to_swlist]->name << "SA" << f->fault_type << " tfsf: " << f->tfsf << " tpsf: " << f->tpsf << " score: " << f->score<< endl;
        }
        return;
    }// if failLog only

    /* test generation mode */
    /* Figure 5 in the PODEM paper */
    while (fault_under_test != nullptr) {
        switch (podem(fault_under_test, current_backtracks)) {
            case TRUE:
                /* form a vector */
                vec.clear();
                for (wptr w: cktin) {
                    vec.push_back(itoc(w->value));
                }
                /*by defect, we want only one pattern per fault */
                /*run a fault simulation, drop ALL detected faults */
                if (total_attempt_num == 1) {
                    fault_sim_a_vector(vec, current_detect_num);
                    total_detect_num += current_detect_num;
                }
                    /* If we want mutiple petterns per fault,
                     * NO fault simulation.  drop ONLY the fault under test */
                else {
                    fault_under_test->detect = TRUE;
                    /* drop fault_under_test */
                    flist_undetect.remove(fault_under_test);
                }
                in_vector_no++;
                break;
            case FALSE:fault_under_test->detect = REDUNDANT;
                no_of_redundant_faults++;
                break;

            case MAYBE:no_of_aborted_faults++;
                break;
        }
        fault_under_test->test_tried = true;
        fault_under_test = nullptr;
        for (fptr fptr_ele: flist_undetect) {
            if (!fptr_ele->test_tried) {
                fault_under_test = fptr_ele;
                break;
            }
        }
        total_no_of_backtracks += current_backtracks; // accumulate number of backtracks
        no_of_calls++;
    }

    display_undetect();
    fprintf(stdout, "\n");
    fprintf(stdout, "#number of aborted faults = %d\n", no_of_aborted_faults);
    fprintf(stdout, "\n");
    fprintf(stdout, "#number of redundant faults = %d\n", no_of_redundant_faults);
    fprintf(stdout, "\n");
    fprintf(stdout, "#number of calling podem1 = %d\n", no_of_calls);
    fprintf(stdout, "\n");
    fprintf(stdout, "#total number of backtracks = %d\n", total_no_of_backtracks);
}/* end of test */


/* constructor of ATPG */
ATPG::ATPG() {
    /* orginally assigned in tpgmain.c */
    this->backtrack_limit = 50;     /* default value */
    this->total_attempt_num = 1;    /* default value */
    this->fsim_only = false;        /* flag to indicate fault simulation only */
    this->tdfsim_only = false;      /* flag to indicate tdfault simulation only */
    this->failLog_only = false;
    this->diag_only = false;
    /* orginally assigned in input.c */
    this->debug = 0;                /* != 0 if debugging;  this is a switch of debug mode */
    this->lineno = 0;               /* current line number */
    this->targc = 0;                /* number of args on current command line */
    this->file_no = 0;              /* number of current file */

    /* orginally assigned in init_flist.c */
    this->num_of_gate_fault = 0; // totle number of faults in the whole circuit

    /* orginally assigned in test.c */
    this->in_vector_no = 0;         /* number of test vectors generated */
}

/* constructor of WIRE */
ATPG::WIRE::WIRE() {
    this->value = 0;
    this->level = 0;
    this->wire_value_g = 0;
    this->wire_value_f = 0;
    this->wlist_index = 0;
}

/* constructor of NODE */
ATPG::NODE::NODE() {
    this->type = 0;
    this->marked = false;
}

/* constructor of FAULT */
ATPG::FAULT::FAULT() {
    this->node = nullptr;
    this->io = 0;
    this->index = 0;
    this->fault_type = 0;
    this->detect = 0;
    this->test_tried = false;
    this->eqv_fault_num = 0;
    this->to_swlist = 0;
    this->fault_no = 0;
    this->tfsf = 0;
    this->tpsf = 0;
    this->tfsp = 0;
    this->score = 0;
}

ATPG::TEST_RESULT::TEST_RESULT() {
    this->node = nullptr;
    this->io = 0;
    this->index = 0;
    this->observed = 0;
    this->expected = 0;
    this->vec = "";
    this->detect = FALSE;
}

