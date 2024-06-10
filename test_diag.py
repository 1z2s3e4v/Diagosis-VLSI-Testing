import pandas as pd
import os, sys
import re
import subprocess

# Define config
atpg_bin = './src/atpg'
patterns_dir = './patterns'
circuits_dir = './sample_circuits'
failLog_dir = './failLog'
diag_rpt_dir = './diag_rpt'
Info_filename = failLog_dir+'/Info_failLog.csv'
circuit_names = ['c17', 'c432', 'c499', 'c880', 'c1355', 'c2670', 'c3540', 'c6288', 'c7552']
fault_names = ['Fault 1', 'Fault 2', 'Fault 3', 'Fault 4', 'Fault 5']
os.system(f"mkdir -p {diag_rpt_dir}")

# Read info_failLog.csv
df = pd.read_csv(Info_filename)
# print(df)

def run_diag(circuit_name, failLog_idx):
    # run atpg (ex: ./src/atpg -diag ./patterns/golden_c17.ptn ./sample_circuits/c17.ckt ./failLog/c17-020.failLog)
    rpt_filename = f'{circuit_name}-{failLog_idx}.log'
    cmd = f"{atpg_bin} -diag {patterns_dir}/golden_{circuit_name}.ptn {circuits_dir}/{circuit_name}.ckt {failLog_dir}/{circuit_name}-{failLog_idx:03}.failLog > {diag_rpt_dir}/{rpt_filename}"
    # print(f"Run cmd: {cmd}")
    time_cmd = f"/usr/bin/time -p {cmd}"
    result = subprocess.run(time_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    # GET time
    user_time = None
    for line in result.stderr.split('\n'):
        if line.startswith("user"):
            user_time = float(re.findall(r'\d+\.\d+', line)[0])
            break
    if user_time is not None:
        return user_time
    else: return -1

def find_fault(golden_f, diag_fs):
    for i in range(len(diag_fs)):
        if golden_f == diag_fs[i]['fault_name'] or golden_f in diag_fs[i]['equivalent_faults']:
            return i
    return -1


def check_diag_rpt(circuit_name, failLog_idx):
    # get diag faults
    rpt_filename = f'{diag_rpt_dir}/{circuit_name}-{failLog_idx}.log'
    with open(rpt_filename, 'r') as file: diag_rpt = file.read()
    extract_pattern = re.compile(r'No\.(\d+)\s+([\w\d()_*]+\s+\w+\s+\w+\s+\w+),\s+groupID=(\d+),\s+TFSF=(\d+),\s+TPSF=(\d+),\s+TFSP=(\d+),\s+TPSP=(\d+),\s+score=([\d.]+)\s+\[ equivalent faults: (.*?)\]')
    matches = extract_pattern.findall(diag_rpt)
    diag_faults = []
    for match in matches:
        diag_f = {}
        # diag_f['rank'] = match[0]
        diag_f['fault_name'] = match[1]; diag_f['group_id'] = match[2]; diag_f['tfsf'] = match[3]; diag_f['tpsf'] = match[4]; diag_f['tfsp'] = match[5]; diag_f['tpsp'] = match[6]; diag_f['score'] = match[7]; diag_f['equivalent_faults'] = [fault.strip() for fault in match[8].strip().split(',') if fault.strip()]
        diag_faults.append(diag_f)
    # get the row in df
    condition = (df['Circuit'] == circuit_name) & (df['Fail Log Index'] == failLog_idx)
    row = df.loc[df[condition].index[0]]
    # get golden faults
    golden_faults = []
    for f in fault_names:
        if row[f] != "None": golden_faults.append(row[f].replace('/', ' '))
    # calculate accuracy and resolution
    num_correctly_diag_faults = 0
    num_total_injected_faults = len(golden_faults)
    num_diag_faults = len(diag_faults)
    diag_accuracy = 0.0
    diag_resolution = 30.0 
    golden_fault_rank_in_diag = []
    for gf in golden_faults:
        diag_rank = find_fault(gf, diag_faults)+1
        if diag_rank <= 5 and diag_rank > 0:
            num_correctly_diag_faults += 1
        golden_fault_rank_in_diag.append(diag_rank)
    diag_accuracy = float(num_correctly_diag_faults) / num_total_injected_faults
    if num_correctly_diag_faults > 0:
        diag_resolution = float(num_diag_faults) / num_correctly_diag_faults
    return diag_accuracy, diag_resolution, [num_correctly_diag_faults, num_total_injected_faults, num_diag_faults, golden_fault_rank_in_diag]

# Check mode (run_one or run_all)
if (len(sys.argv) == 2 and sys.argv[1] == "all"): # Run all
    # Check for each row
    for c_name in circuit_names:
        diag_accs = []
        diag_resols = []
        diag_runtimes = []
        # Check for each row
        for index, row in df.iterrows():
            # get info
            circuit_name = row["Circuit"]
            failLog_idx = row["Fail Log Index"]
            # Check for round
            if circuit_name != c_name: continue
            # run atpg 
            runtime = run_diag(circuit_name, failLog_idx)
            # check diagnosis accuracy
            diag_accuracy, diag_resolution, other_info = check_diag_rpt(circuit_name, failLog_idx)
            diag_accs.append(diag_accuracy); diag_resols.append(diag_resolution); diag_runtimes.append(runtime); 
            print(f"{circuit_name}-{failLog_idx:03}: Acc={diag_accuracy:.2f} ({other_info[0]}/{other_info[1]}), Resol={diag_resolution:.2f} ({other_info[2]}/{other_info[0]}), Runtime={runtime:.2f}s", end='')
            print(", diag_ranks of injected faults: ", end='')
            for rank in other_info[3]: print(f" {rank},", end='')
            print("")
        print("------------------------------------------------------------------")
        print(f"{c_name} summary: Avg Acc={sum(diag_accs)/len(diag_accs):.2f}, Avg Resol={sum(diag_resols)/len(diag_resols):.2f}, Avg Runtime={sum(diag_runtimes)/len(diag_runtimes):.2f}s")
        print("")

elif len(sys.argv) == 3: # Run one case
    # run atpg 
    runtime = run_diag(sys.argv[1], int(sys.argv[2]))
    # check diagnosis accuracy
    diag_accuracy, diag_resolution, other_info = check_diag_rpt(sys.argv[1], int(sys.argv[2]))
    print(f"Acc={diag_accuracy:.2f} ({other_info[0]}/{other_info[1]}), Resol={diag_resolution:.2f} ({other_info[2]}/{other_info[0]}), Runtime={runtime:.2f}s", end='')
    print(", diag_ranks of injected faults: ", end='')
    for rank in other_info[3]: print(f" {rank},", end='')
    print("")

elif len(sys.argv) == 2 and sys.argv[1] in circuit_names: # Run all failLogs for one circuit
    diag_accs = []
    diag_resols = []
    diag_runtimes = []
    # Check for each row
    for index, row in df.iterrows():
        # get info
        circuit_name = row["Circuit"]
        if circuit_name != sys.argv[1]: continue
        failLog_idx = row["Fail Log Index"]
        # run atpg 
        runtime = run_diag(circuit_name, failLog_idx)
        # check diagnosis accuracy
        diag_accuracy, diag_resolution, other_info = check_diag_rpt(circuit_name, failLog_idx)
        diag_accs.append(diag_accuracy); diag_resols.append(diag_resolution); diag_runtimes.append(runtime); 
        print(f"{circuit_name}-{failLog_idx:03}: Acc={diag_accuracy:.2f} ({other_info[0]}/{other_info[1]}), Resol={diag_resolution:.2f} ({other_info[2]}/{other_info[0]}), Runtime={runtime:.2f}s", end='')
        print(", diag_ranks of injected faults: ", end='')
        for rank in other_info[3]: print(f" {rank},", end='')
        print("")
    print("------------------------------------------------------------------")
    print(f"Summary: Avg Acc={sum(diag_accs)/len(diag_accs):.2f}, Avg Resol={sum(diag_resols)/len(diag_resols):.2f}, Avg Runtime={sum(diag_runtimes)/len(diag_runtimes):.2f}s")
    print("")

else:
    print("Usage1: python test_diag.py all")
    print("Usage2: python test_diag.py <circuit_name> <failLog_idx>  (ex: python test_diag.py c17 20)")
    print("Usage2: python test_diag.py <circuit_name>                (ex: python test_diag.py c17)")

