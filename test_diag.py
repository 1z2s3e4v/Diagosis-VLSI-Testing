import pandas as pd
import os, sys
import re
import subprocess
import csv

# Define config
atpg_bin = './src/atpg'
#atpg_bin = './bin_reference/atpg_reference'
patterns_dir = './patterns'
circuits_dir = './sample_circuits'
failLog_dir = './failLog'
diag_rpt_dir = './diag_rpt'
Info_filename = failLog_dir+'/Info_failLog.csv'
circuit_names = ['c17', 'c432', 'c499', 'c880', 'c1355', 'c2670', 'c3540', 'c6288', 'c7552']
# circuit_names = ['c17', 'c432', 'c499', 'c880', 'c1355', 'c2670', 'c3540', 'c7552']
fault_names = ['Fault 1', 'Fault 2', 'Fault 3', 'Fault 4', 'Fault 5']
os.system(f"mkdir -p {diag_rpt_dir}")

# Read info_failLog.csv
df = pd.read_csv(Info_filename)
# print(df)

def run_genFailLog(circuit_name, failLog_idx, diag_faults):
    # run atpg (ex: ./src/atpg -genFailLog ./patterns/golden_c499.ptn ./sample_circuits/c499.ckt -fault ID7"("7")" g389 GI SA1 -fault ID16"("16")" g52 GI SA1)
    failLog_filename = f'{diag_rpt_dir}/{circuit_name}-{failLog_idx}.failLog'
    cmd = f"{atpg_bin} -genFailLog {patterns_dir}/golden_{circuit_name}.ptn {circuits_dir}/{circuit_name}.ckt"
    for diag_f in diag_faults:
        diag_f_info = diag_f.replace('(','\"(\"').replace(')','\")\"')
        cmd += f" -fault {diag_f_info}"
    cmd += f" > {failLog_filename}"

    # print(f"Run cmd: {cmd}")
    subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)

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

def check_failLog(circuit_name, failLog_idx):
    # get failLogs
    golden_failLog_filename = f'{failLog_dir}/{circuit_name}-{failLog_idx:03}.failLog'
    diag_failLog_filename = f'{diag_rpt_dir}/{circuit_name}-{failLog_idx}.failLog'
    golden_failLogs = []
    diag_failLogs = []
    with open(golden_failLog_filename, 'r') as file:
        golden_failLogs = [line.strip() for line in file if line.startswith("vector")]
    with open(diag_failLog_filename, 'r') as file:
        diag_failLogs = [line.strip() for line in file if line.startswith("vector")]
    # check if golden_failLog same with diag_failLog
    missing_vectors = [fail for fail in golden_failLogs if fail not in diag_failLogs]
    perfect_matched = True
    if len(golden_failLogs) != len(diag_failLogs): perfect_matched = False
    for fail_idx, miss in enumerate(missing_vectors):
        if miss: perfect_matched = False#; print(f"  Missed failLog: {golden_failLogs[fail_idx]}");
    return perfect_matched

def check_diag_rpt(circuit_name, failLog_idx):
    # get diag faults
    rpt_filename = f'{diag_rpt_dir}/{circuit_name}-{failLog_idx}.log'
    with open(rpt_filename, 'r') as file: diag_rpt = file.read()
    extract_pattern = re.compile(r'No\.(\d+)\s+([\w\d()_*]+\s+\w+\s+\w+\s+\w+),\s+(.*)\s+\[ equivalent faults: (.*?)\]')
    matches = extract_pattern.findall(diag_rpt)
    diag_faults = []
    for match in matches:
        diag_f = {}
        # diag_f['rank'] = match[0]
        diag_f['fault_name'] = match[1]; diag_f['equivalent_faults'] = [fault.strip() for fault in match[3].strip().split(',') if fault.strip()]
        diag_faults.append(diag_f)
    # get the row in df
    condition = (df['Circuit'] == circuit_name) & (df['Fail Log Index'] == failLog_idx)
    row = df.loc[df[condition].index[0]]
    # get golden faults
    golden_faults = []
    for f in fault_names:
        f_name = str(row[f])
        if f_name != "None" and f_name != "NaN" and f_name != "none" and f_name != "nan" and f_name != "Nan": 
            golden_faults.append(str(row[f]).replace('/', ' '))
    # calculate accuracy and resolution
    num_correctly_diag_faults = 0
    num_total_injected_faults = len(golden_faults)
    num_diag_faults = len(diag_faults)
    diag_accuracy = 0.0
    diag_resolution = 30.0 
    golden_fault_rank_in_diag = []
    # 1. Check the failLog first
    perfect_matched = False
    run_genFailLog(circuit_name, failLog_idx, [diag_f['fault_name'] for diag_f in diag_faults])
    perfect_matched = check_failLog(circuit_name, failLog_idx)
    if perfect_matched:
        diag_accuracy = 1.0
        diag_resolution = 1.0
        num_correctly_diag_faults = len(diag_faults)
        num_diag_faults = len(diag_faults)
    # 2. Check the diag faults if the failLog is not same
    else: 
        for gf in golden_faults:
            diag_rank = find_fault(gf, diag_faults)+1
            if diag_rank <= 5 and diag_rank > 0:
                num_correctly_diag_faults += 1
            golden_fault_rank_in_diag.append(diag_rank)
        diag_accuracy = float(num_correctly_diag_faults) / num_total_injected_faults
        if num_correctly_diag_faults > 0:
            diag_resolution = float(num_diag_faults) / num_correctly_diag_faults
    # return
    other_info = {"num_correctly_diag_faults":num_correctly_diag_faults, "num_total_injected_faults":num_total_injected_faults, "num_diag_faults":num_diag_faults, "golden_fault_rank_in_diag":golden_fault_rank_in_diag, "perfect_matched":perfect_matched}
    return diag_accuracy, diag_resolution, other_info

# Check mode (run_one or run_all)
if (len(sys.argv) == 2 and sys.argv[1] == "all"): # Run all
    csv_columns = ["Circuit Name", "Fail Log Index", "Accuracy", "Accuracy Info", "Resolution", "Resolution Info", "Runtime", "Injected Fault Ranks"]
    with open('results.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(csv_columns)
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
            print(f"{circuit_name}-{failLog_idx:03}: Acc={diag_accuracy:.2f} ({other_info['num_correctly_diag_faults']}/{other_info['num_total_injected_faults']}), Resol={diag_resolution:.2f} ({other_info['num_diag_faults']}/{other_info['num_correctly_diag_faults']}), Runtime={runtime:.2f}s", end='')
            if other_info["perfect_matched"]:
                print(", failLog perfectly matched!")
            else:
                print(", diag_ranks of injected faults: ", end='')
                for rank in other_info["golden_fault_rank_in_diag"]: print(f" {rank},", end='')
                print("")
            # output csv
            accuracy_info = f"'{other_info['num_correctly_diag_faults']}/{other_info['num_total_injected_faults']}"
            resolution_info = f"'{other_info['num_diag_faults']}/{other_info['num_correctly_diag_faults']}"
            injected_fault_ranks = ", ".join(map(str, other_info["golden_fault_rank_in_diag"]))
            row_data = [circuit_name, f"{failLog_idx}", f"{diag_accuracy}", accuracy_info, f"{diag_resolution}", resolution_info, f"{runtime}", injected_fault_ranks]
            with open('results.csv', mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row_data)
        print("------------------------------------------------------------------")
        print(f"{c_name} summary: Avg Acc={sum(diag_accs)/len(diag_accs):.2f}, Avg Resol={sum(diag_resols)/len(diag_resols):.2f}, Avg Runtime={sum(diag_runtimes)/len(diag_runtimes):.2f}s")
        print("")

elif len(sys.argv) == 3: # Run one case
    # run atpg 
    runtime = run_diag(sys.argv[1], int(sys.argv[2]))
    # check diagnosis accuracy
    diag_accuracy, diag_resolution, other_info = check_diag_rpt(sys.argv[1], int(sys.argv[2]))
    print(f"Acc={diag_accuracy:.2f} ({other_info['num_correctly_diag_faults']}/{other_info['num_total_injected_faults']}), Resol={diag_resolution:.2f} ({other_info['num_diag_faults']}/{other_info['num_correctly_diag_faults']}), Runtime={runtime:.2f}s", end='')
    if other_info["perfect_matched"]:
        print(", failLog perfectly matched!")
    else:
        print(", diag_ranks of injected faults: ", end='')
        for rank in other_info["golden_fault_rank_in_diag"]: print(f" {rank},", end='')
        print("")

elif len(sys.argv) == 2 and sys.argv[1] in circuit_names: # Run all failLogs for one circuit
    diag_accs = []
    diag_resols = []
    diag_runtimes = []
    csv_columns = ["Circuit Name", "Fail Log Index", "Accuracy", "Accuracy Info", "Resolution", "Resolution Info", "Runtime", "Injected Fault Ranks"]
    with open(f'results-{sys.argv[1]}.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(csv_columns)
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
        print(f"{circuit_name}-{failLog_idx:03}: Acc={diag_accuracy:.2f} ({other_info['num_correctly_diag_faults']}/{other_info['num_total_injected_faults']}), Resol={diag_resolution:.2f} ({other_info['num_diag_faults']}/{other_info['num_correctly_diag_faults']}), Runtime={runtime:.2f}s", end='')
        if other_info["perfect_matched"]:
            print(", failLog perfectly matched!")
        else:
            print(", diag_ranks of injected faults: ", end='')
            for rank in other_info["golden_fault_rank_in_diag"]: print(f" {rank},", end='')
            print("")
        # output csv
        accuracy_info = f"'{other_info['num_correctly_diag_faults']}/{other_info['num_total_injected_faults']}"
        resolution_info = f"'{other_info['num_diag_faults']}/{other_info['num_correctly_diag_faults']}"
        injected_fault_ranks = ", ".join(map(str, other_info["golden_fault_rank_in_diag"]))
        row_data = [circuit_name, f"{failLog_idx}", f"{diag_accuracy}", accuracy_info, f"{diag_resolution}", resolution_info, f"{runtime}", injected_fault_ranks]
        with open('results.csv', mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row_data)
    print("------------------------------------------------------------------")
    print(f"Summary: Avg Acc={sum(diag_accs)/len(diag_accs):.2f}, Avg Resol={sum(diag_resols)/len(diag_resols):.2f}, Avg Runtime={sum(diag_runtimes)/len(diag_runtimes):.2f}s")
    print("")

else:
    print("Usage1: python test_diag.py all")
    print("Usage2: python test_diag.py <circuit_name> <failLog_idx>  (ex: python test_diag.py c17 20)")
    print("Usage2: python test_diag.py <circuit_name>                (ex: python test_diag.py c17)")

