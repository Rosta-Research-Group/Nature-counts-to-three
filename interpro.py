import os
import argparse
import re
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 16})


def get_ipr_per_chain():
    ipr = glob("data/ipr/*dat")
    supfam = ["IPR000336", "IPR000415", "IPR003789", "IPR004115", "IPR008210", "IPR008250", "IPR008916", "IPR008919",
              "IPR008921", "IPR008925", "IPR008936", "IPR008948", "IPR008978", "IPR008979", "IPR008984", "IPR008988",
              "IPR008995", "IPR009000", "IPR009001", "IPR009003", "IPR009008", "IPR009010", "IPR009030", "IPR009050",
              "IPR009057", "IPR009060", "IPR009061", "IPR009068", "IPR009078", "IPR009080", "IPR009091", "IPR009097",
              "IPR010978", "IPR010982", "IPR010994", "IPR010997", "IPR010999", "IPR011004", "IPR011005", "IPR011006",
              "IPR011008", "IPR011009", "IPR011011", "IPR011026", "IPR011029", "IPR011035", "IPR011037", "IPR011042",
              "IPR011043", "IPR011044", "IPR011047", "IPR011049", "IPR011050", "IPR011051", "IPR011053", "IPR011054",
              "IPR011059", "IPR011060", "IPR011068", "IPR011162", "IPR011254", "IPR011257", "IPR011330", "IPR011335",
              "IPR011989", "IPR011990", "IPR011992", "IPR011993", "IPR012292", "IPR012334", "IPR012337", "IPR012340",
              "IPR012344", "IPR012348", "IPR012675", "IPR012676", "IPR012677", "IPR013035", "IPR013083", "IPR013320",
              "IPR013755", "IPR013756", "IPR013761", "IPR013783", "IPR013784", "IPR013785", "IPR013792", "IPR013806",
              "IPR013815", "IPR014042", "IPR014049", "IPR014352", "IPR014709", "IPR014710", "IPR014721", "IPR014724",
              "IPR014729", "IPR014746", "IPR014756", "IPR015421", "IPR015422", "IPR015424", "IPR015797", "IPR015806",
              "IPR015813", "IPR015824", "IPR015915", "IPR015919", "IPR015943", "IPR015947", "IPR015955", "IPR016024",
              "IPR016039", "IPR016064", "IPR016084", "IPR016102", "IPR016120", "IPR016135", "IPR016142", "IPR016143",
              "IPR016155", "IPR016161", "IPR016162", "IPR016163", "IPR016181", "IPR016185", "IPR016186", "IPR016187",
              "IPR016193", "IPR016195", "IPR017437", "IPR017438", "IPR017449", "IPR017850", "IPR017856", "IPR017945",
              "IPR017946", "IPR018162", "IPR018163", "IPR018193", "IPR018197", "IPR018490", "IPR020056", "IPR020061",
              "IPR020568", "IPR020751", "IPR020752", "IPR020825", "IPR021109", "IPR022384", "IPR022636", "IPR022830",
              "IPR023098", "IPR023142", "IPR023168", "IPR023183", "IPR023191", "IPR023198", "IPR023211", "IPR023213",
              "IPR023214", "IPR023293", "IPR023298", "IPR023299", "IPR023314", "IPR023318", "IPR023366", "IPR023382",
              "IPR023420", "IPR023465", "IPR023509", "IPR023578", "IPR023580", "IPR024034", "IPR024073", "IPR024074",
              "IPR024075", "IPR024079", "IPR024096", "IPR024185", "IPR027267", "IPR027409", "IPR027410", "IPR027413",
              "IPR027417", "IPR027421", "IPR027434", "IPR027442", "IPR027460", "IPR027483", "IPR027484", "IPR028082",
              "IPR028375", "IPR029000", "IPR029001", "IPR029016", "IPR029021", "IPR029033", "IPR029038", "IPR029040",
              "IPR029044", "IPR029045", "IPR029047", "IPR029048", "IPR029052", "IPR029053", "IPR029055", "IPR029056",
              "IPR029057", "IPR029058", "IPR029060", "IPR029062", "IPR029063", "IPR029067", "IPR029068", "IPR029071",
              "IPR029151", "IPR029787", "IPR032466", "IPR032675", "IPR032710", "IPR033469", "IPR035073", "IPR035099",
              "IPR035437", "IPR035891", "IPR035892", "IPR035896", "IPR035899", "IPR035902", "IPR035907", "IPR035911",
              "IPR035913", "IPR035919", "IPR035921", "IPR035929", "IPR035930", "IPR035959", "IPR035963", "IPR035965",
              "IPR035966", "IPR035979", "IPR035983", "IPR035985", "IPR035990", "IPR036020", "IPR036024", "IPR036025",
              "IPR036028", "IPR036034", "IPR036043", "IPR036046", "IPR036047", "IPR036055", "IPR036061", "IPR036075",
              "IPR036086", "IPR036097", "IPR036116", "IPR036117", "IPR036121", "IPR036127", "IPR036129", "IPR036135",
              "IPR036137", "IPR036144", "IPR036157", "IPR036161", "IPR036163", "IPR036179", "IPR036181", "IPR036193",
              "IPR036203", "IPR036206", "IPR036213", "IPR036225", "IPR036236", "IPR036237", "IPR036249", "IPR036253",
              "IPR036258", "IPR036264", "IPR036265", "IPR036266", "IPR036274", "IPR036279", "IPR036282", "IPR036291",
              "IPR036305", "IPR036320", "IPR036322", "IPR036332", "IPR036333", "IPR036352", "IPR036371", "IPR036380",
              "IPR036388", "IPR036390", "IPR036393", "IPR036397", "IPR036409", "IPR036412", "IPR036420", "IPR036425",
              "IPR036426", "IPR036427", "IPR036430", "IPR036451", "IPR036457", "IPR036465", "IPR036477", "IPR036480",
              "IPR036481", "IPR036499", "IPR036513", "IPR036514", "IPR036522", "IPR036523", "IPR036526", "IPR036551",
              "IPR036553", "IPR036554", "IPR036565", "IPR036566", "IPR036571", "IPR036572", "IPR036597", "IPR036599",
              "IPR036603", "IPR036604", "IPR036612", "IPR036615", "IPR036621", "IPR036637", "IPR036640", "IPR036641",
              "IPR036643", "IPR036644", "IPR036651", "IPR036654", "IPR036670", "IPR036674", "IPR036676", "IPR036688",
              "IPR036690", "IPR036691", "IPR036695", "IPR036736", "IPR036738", "IPR036741", "IPR036745", "IPR036754",
              "IPR036759", "IPR036768", "IPR036770", "IPR036775", "IPR036779", "IPR036784", "IPR036790", "IPR036802",
              "IPR036830", "IPR036832", "IPR036844", "IPR036850", "IPR036855", "IPR036860", "IPR036862", "IPR036865",
              "IPR036866", "IPR036867", "IPR036869", "IPR036871", "IPR036872", "IPR036873", "IPR036875", "IPR036888",
              "IPR036890", "IPR036891", "IPR036892", "IPR036895", "IPR036897", "IPR036898", "IPR036901", "IPR036907",
              "IPR036913", "IPR036914", "IPR036915", "IPR036918", "IPR036921", "IPR036925", "IPR036928", "IPR036936",
              "IPR036940", "IPR036941", "IPR036945", "IPR036946", "IPR036947", "IPR036957", "IPR036961", "IPR036964",
              "IPR036968", "IPR036969", "IPR036974", "IPR036986", "IPR036988", "IPR037006", "IPR037009", "IPR037013",
              "IPR037017", "IPR037033", "IPR037034", "IPR037035", "IPR037051", "IPR037052", "IPR037064", "IPR037070",
              "IPR037080", "IPR037085", "IPR037118", "IPR037123", "IPR037136", "IPR037151", "IPR037159", "IPR037160",
              "IPR037162", "IPR037171", "IPR037172", "IPR037176", "IPR037196", "IPR037204", "IPR037214", "IPR037227",
              "IPR037230", "IPR037238", "IPR037243", "IPR037257", "IPR037265", "IPR038030", "IPR038055", "IPR038083",
              "IPR038120", "IPR038123", "IPR038124", "IPR038138", "IPR038152", "IPR038154", "IPR038155", "IPR038158",
              "IPR038166", "IPR038170", "IPR038178", "IPR038188", "IPR038199", "IPR038222", "IPR038232", "IPR038238",
              "IPR038239", "IPR038249", "IPR038252", "IPR038256", "IPR038286", "IPR038302", "IPR038318", "IPR038320",
              "IPR038321", "IPR038331", "IPR038337", "IPR038345", "IPR038346", "IPR038357", "IPR038376", "IPR038385",
              "IPR038388", "IPR038400", "IPR038408", "IPR038411", "IPR038418", "IPR038419", "IPR038424", "IPR038428",
              "IPR038429", "IPR038451", "IPR038469", "IPR038515", "IPR038557", "IPR038568", "IPR038593", "IPR038607",
              "IPR038614", "IPR038634", "IPR038662", "IPR038677", "IPR038688", "IPR038718", "IPR038763", "IPR038765",
              "IPR040442", "IPR041856", "IPR041931", "IPR041948", "IPR042032", "IPR042035", "IPR042063", "IPR042078",
              "IPR042079", "IPR042087", "IPR042090", "IPR042095", "IPR042099", "IPR042101", "IPR042102", "IPR042103",
              "IPR042107", "IPR042109", "IPR042110", "IPR042111", "IPR042176", "IPR042199", "IPR042205", "IPR042209",
              "IPR042213", "IPR042236", "IPR042240", "IPR042257", "IPR042258", "IPR042267", "IPR042295", "IPR042302",
              "IPR042308", "IPR042309", "IPR042310", "IPR042311", "IPR042449", "IPR042463", "IPR042515", "IPR042521",
              "IPR042542", "IPR042543", "IPR042544", "IPR042558", "IPR042559", "IPR042570", "IPR043024", "IPR043056",
              "IPR043110", "IPR043119", "IPR043120", "IPR043128", "IPR043129", "IPR043130", "IPR043133", "IPR043134",
              "IPR043136", "IPR043150", "IPR043171", "IPR043174", "IPR043177", "IPR043178", "IPR043181", "IPR043472",
              "IPR043477", "IPR043478", "IPR043502", "IPR043503", "IPR043504", "IPR043519", "IPR044876", "IPR044893",
              "IPR044894", "IPR044896", "IPR044912", "IPR044923", "IPR044925", "IPR044926", "IPR044929", "IPR044936",
              "IPR044949", "IPR045851", "IPR045860", "IPR045864", "IPR045865", "IPR046342", "IPR046346", "IPR046348",
              "IPR046349", "IPR046357", "IPR046375", "IPR046409", "IPR046428", "IPR046429", "IPR046430", "IPR046437",
              "IPR046812", "IPR046813", "IPR046814", "IPR046933", "IPR046938", "IPR047108"]
    ipr_per_chain = {}
    for i in ipr:
        id = i.split('_')[0]
        with open(i, "r") as f:
            for line in f:
                l = line.strip().split(' ')
                if len(l) != 3:
                    chains = ""
                    for x in line[4:].strip():
                        try:
                            int(x)
                            break
                        except ValueError:
                            chains += x
                    l = [line[:4], chains]
                for c in l[1].split(','):
                    if len(c) == 0:
                        continue
                    if f"{l[0]}_{c}" not in ipr_per_chain:
                        ipr_per_chain[f"{l[0]}_{c}"] = {"sf": [], "f": []}
                    if id in supfam:
                        ipr_per_chain[f"{l[0]}_{c}"]["sf"].append(id)
                    else:
                        ipr_per_chain[f"{l[0]}_{c}"]["f"].append(id)
    with open("data/fams.csv", "w") as fout:
        for c in ipr_per_chain:
            fout.write(f"{c},{len(ipr_per_chain[c]['f'])}")
            for f in ipr_per_chain[c]["f"]:
                fout.write(f",{f}")
            fout.write("\n")
    with open("data/homSF.csv", "w") as fout:
        for c in ipr_per_chain:
            fout.write(f"{c},{len(ipr_per_chain[c]['sf'])}")
            for f in ipr_per_chain[c]["sf"]:
                fout.write(f",{f}")
            fout.write("\n")
    return ipr_per_chain


def get_ipr_corr(ipr_per_chain):
    fam = []
    hsf = []
    for k in ipr_per_chain:
        fam += ipr_per_chain[k]["f"]
        hsf += ipr_per_chain[k]["sf"]
    fam = np.unique(fam)
    hsf = np.unique(hsf)
    fammap = np.zeros(shape=(len(ipr_per_chain), len(fam)), dtype=bool)
    hsfmap = np.zeros(shape=(len(ipr_per_chain), len(hsf)), dtype=bool)
    i = 0
    for k in ipr_per_chain:
        for j in range(len(ipr_per_chain[k]["f"])):
            fammap[i, np.where(fam == ipr_per_chain[k]["f"][j])] = True
        for j in range(len(ipr_per_chain[k]["sf"])):
            hsfmap[i, np.where(hsf == ipr_per_chain[k]["sf"][j])] = True
        i += 1
    famcorr = np.corrcoef(fammap, rowvar=False)
    hsfcorr = np.corrcoef(hsfmap, rowvar=False)
    xcorr = np.corrcoef(fammap, hsfmap, rowvar=False)[:len(fam), len(fam):]
    sns.heatmap(famcorr)
    plt.tight_layout()
    plt.savefig("famcorr.png", dpi=600)
    plt.close()
    sns.heatmap(hsfcorr)
    plt.tight_layout()
    plt.savefig("hsfcorr.png", dpi=600)
    plt.close()
    sns.heatmap(xcorr)
    plt.tight_layout()
    plt.savefig("xcorr.png", dpi=600)
    plt.close()
    return


def parse_ipr_ec(datafile="interpro.xml"):
    """
    The datafile should be downloaded from the interpro website.
    :param datafile:
    :type datafile:
    :return:
    :rtype:
    """
    ec = {}
    with open(datafile, "r") as f:
        doclist = False
        cur_id = ""
        for line in f:
            if "<interpro id=" in line:
                cur_id = line.split('"')[1].split(">")[0]
                if cur_id not in ec:
                    ec[cur_id] = []
            elif "<external_doc_list>" in line:
                doclist = True
            elif  "</external_doc_list>" in line:
                doclist = False
            elif doclist and 'db="EC"' in line:
                this_ec = line.split('"')[3]
                if '-' not in this_ec and this_ec not in ec[cur_id]:
                    ec[cur_id].append(this_ec)
    return ec


def get_chain_idx(pdbfile):
    code = {}
    i = 0
    with open(pdbfile, "r") as f:
        for line in f:
            if "TER" == line[:3]:
                code[line[21]] = i
                i += 1
    return code


def process_ipr(iprnum, pdb_chain_supfam, output=None):
    fout = open(output, "w")
    with (open(f"data/ipr/{iprnum}_chain.dat", "r") as f):
        for line in f:
            l = line.strip().split(' ')
            # if len(pdb_chain_supfam[pdb_chain_supfam.PDB == l[0].upper()]) == 0:
            #     print(f"skipping {l[0]}, did not go to supfam")
            #     continue
            chain_ids = l[1].split(',')
            try:
                ranges = l[2].split(',')
                for i in range(len(ranges)):
                    if ranges[i] == '':
                        ranges[i] = "0..9999"
            except IndexError:
                ranges = ["0..9999" for _ in chain_ids]
            if len(chain_ids) < len(ranges):
                new_ranges = []
                old = np.inf
                for x in ranges:
                    lims = x.split('.')
                    if int(lims[0]) < old:
                        new_ranges.append(f"{x}")
                    else:
                        new_ranges[-1] += f"-{x}"
                    old = int(lims[-1])
                # print("variable number of ranges")
                if len(new_ranges) != len(chain_ids):
                    ranges = ["0..9999" for _ in chain_ids]
                else:
                    ranges = new_ranges
            for i in range(len(chain_ids)):
                c = chain_ids[i]
                if len(c) == 0:
                    continue
                view = pdb_chain_supfam[pdb_chain_supfam.PDB == l[0].upper()]
                if len(view[view['chain_ipro'] == c]) == 0:
                    if output is None:
                        print(f"{l[0]}_{c} is not found in supfam data")
                    else:
                        fout.write(f"{l[0]}_{c} is not found in supfam data\n")
                    continue
                else:
                    resp = f"{l[0]}_{c};"
                    resp += f"{ranges[i]}:"
                    resp += "supfam hits;range:"
                    view = view[view['chain_ipro'] == c]
                    for j in view.index:
                        resp += "{0};{1}..{2}:".format(view.Superfamily[j],
                                                      view["Start_Residue"][j], view["End_Residue"][j])
                    if output is None:
                        print(resp)
                    else:
                        fout.write(resp)
                        fout.write("\n")
    fout.close()
    return


def count_chains(iprnum):
    count = 0
    with open(f"{iprnum}_chain.dat", "r") as f:
        for line in f:
            l = line.strip().split(' ')
            for c in l[1].split(','):
                if len(c) != 0:
                    count += 1
    print(f"{iprnum}: {count}")
    return


def is_part_of(sfr, iprr):
    r1 = [[int(x.split('.')[0]), int(x.split('.')[-1])] for x in sfr.split('-')]
    c1 = [np.mean(x) for x in r1]
    # s1 = [x[1] - x[0] for x in r1]
    r2 = [[int(x.split('.')[0]), int(x.split('.')[-1])] for x in iprr.split('-')]
    # c2 = [np.mean(x) for x in r2]
    s2 = [x[1] - x[0] for x in r2]
    for i in range(len(r1)):
        for j in range(len(r2)):
            if c1[i] >= r2[j][0] and c1[i] <= r2[j][1]:
                return True
            # if s1[i] >= s2[j] and c2[j] >= r1[i][0] and c2[j] <= r1[i][1]:
            #     return True
            # elif s1[i] <= s2[j] and c1[i] >= r2[j][0] and c1[i] <= r2[j][1]:
            #     return True
    return False


def parse_ipr_proc(iprnum, suffix="_proc.log"):
    with open(f"{iprnum}{suffix}", "r") as f:
        lines = f.readlines()
    pdbc = [line.split(';')[0].split()[0].replace("_", "") for line in lines]
    sites = {"p1": 0, "p2": 0, "p3": 0}
    pdbs = {"p1": [], "p2": [], "p3": [], "no_P": [], "no_data": []}
    for pc in pdbc:
        try:
            p = phosphate["numP"][pc.upper()]
            sites[f"p{p}"] += 1
            if pc not in pdbs[f"p{p}"]:
                pdbs[f"p{p}"].append(pc)
        except KeyError:
            try:
                if pc not in pdbs["no_P"]:
                    pdbs["no_P"].append(pc)
            except KeyError:
                if pc not in pdbs["no_data"]:
                    pdbs["no_data"].append(pc)
    if sites["p1"] + sites["p2"] + sites["p3"] == 0:
        return
    supfams_in = {}
    in_supfams = {}
    out_supfams = {}
    hitchains = []
    for line in lines:
        if "hits" in line:
            chain, ranges = line.split(":")[0].split(';')
            hitchains.append(chain)
            sfr = line.strip().split(":")[2:-1]
            sf = [x.split(';')[0] for x in sfr]
            supfamr = [x.split(';')[1] for x in sfr]
            sfin = []
            iprin = []
            out = []
            for i in range(len(sf)):
                if is_part_of(supfamr[i], ranges):
                    if sf[i] not in sfin:
                        sfin.append(sf[i])
                        try:
                            supfams_in[sf[i]] += 1
                        except KeyError:
                            supfams_in[sf[i]] = 1
                elif is_part_of(ranges, supfamr[i]):
                    if sf[i] not in iprin:
                        iprin.append(sf[i])
                        try:
                            in_supfams[sf[i]] += 1
                        except KeyError:
                            in_supfams[sf[i]] = 1
                elif sf[i] not in out:
                    out.append(sf[i])
                    try:
                        out_supfams[sf[i]] += 1
                    except KeyError:
                        out_supfams[sf[i]] = 1
    repr = 0
    single_repr_count = 0
    repr_entry = ""
    for k in supfams_in:
        if k in repr_supfam:
            repr += 1
            single_repr_count = supfams_in[k]
            repr_entry = k
    if repr == 1 and single_repr_count / len(pdbc) >= 0.5:
        # print(f"{iprnum} {repr_entry}")
        pass
    elif len(hitchains) / len(pdbc) < 0.5:
        pass
    else:
        pass
    print(
        f"https://www.ebi.ac.uk/interpro/entry/InterPro/{iprnum};{len(hitchains)};{len(pdbc)};{sites['p1']};{sites['p2']};{sites['p3']};{' '.join(pdbs['p3'])};{' '.join(pdbs['p2'])};{' '.join(pdbs['p1'])};{' '.join(pdbs['no_P'])};{' '.join(pdbs['no_data'])}")
    for k in supfams_in:
        print(f";{k};{supfams_in[k]};SUPFAM IN IPR")
    for k in in_supfams:
        print(f";{k};{in_supfams[k]};IPR IN SUPFAM")
    for k in out_supfams:
        print(f";{k};{out_supfams[k]};OUT")

    return


def find_chain_singular_IPR_supfam(iprnum, suffix="_proc.log"):
    iprstats = pd.read_csv("data/ipr_chain.csv", sep=';', index_col=0)
    with open(f"{iprnum}{suffix}", "r") as f:
        for line in f:
            if "hits" in line:
                c = line.split(" ")[0]
                # if iprstats["hsf"][c] + iprstats["fam"][c] > 1:
                #     print(f"{c};covers multiple IPR numbers")
                #     continue
                supfams = []
                sf = line.strip().split(":")[1:]
                for s in sf:
                    if len(s) < 2:
                        continue
                    if s not in supfams:
                        supfams.append(s)
                if len(supfams) > 1:
                    print(f"{c};covers multiple supfams")
                    continue
                elif len(supfams) == 1:
                    if iprstats["hsf"][c] == 1:
                        print(f"{c};{iprnum};{supfams[0]};IPR unique superfamily;besides {iprstats['fam'][c]} family")
                    elif iprstats["fam"][c] == 1:
                        print(f"{c};{iprnum};{supfams[0]};IPR unique family;besides {iprstats['hsf'][c]} superfamily")
    return


def get_fasta(iprnum, pdb_chain_supfam):
    if not os.path.isfile("pdb_seqres.txt"):
        print("pdb_seqres.txt should be downloaded from the RCSB website to work with the entire PDB")
        raise FileNotFoundError
    seqs = {}
    with open("pdb_seqres.txt", "r") as f:
        for line in f:
            if line.startswith(">"):
                c = line[1:].split(' ')[0]
            else:
                seqs[c] = line.strip()
    with open(f"{iprnum}_chain.dat", "r") as f:
        fout = open(f"{iprnum}.fa", "w")
        for line in f:
            l = line.strip().split(' ')
            view = pdb_chain_supfam[pdb_chain_supfam.PDB == l[0].upper()]
            chain_ids = l[1].split(',')
            try:
                ranges = l[2].split(',')
            except IndexError:
                ranges = ["0..9999" for _ in chain_ids]
            for c in chain_ids:
                if len(c) == 0:
                    continue
                if len(view[view['chain_ipro'] == c]) > 0:
                    print(f"{l[0]}_{c} is already processed")
                try:
                    fout.write(f">{l[0]}_{c}\n{seqs[f'{l[0]}_{c}']}\n")
                except KeyError:
                    print(f"{l[0]}_{c} is not found in rcsb data")
        fout.close()
    return


def chain_labeling(pdb_chain_supfam):
    if not os.path.isfile("pdb_seqres.txt"):
        print("pdb_seqres.txt should be downloaded from the RCSB website to work with the entire PDB")
        raise FileNotFoundError
    pdb_seq_list = []
    with open("pdb_seqres.txt", "r") as f:
        for line in f:
            if line[0] == ">":
                pdb_seq_list.append(line.strip())
    pdb_chains = {}
    for l in pdb_seq_list:
        p = l[1:5]
        if p not in pdb_chains:
            pdb_chains[p] = []
        pdb_chains[p].append(l.split(" ")[0].split('_')[1])
    for i in range(len(pdb_chain_supfam)):
        if ~np.isnan(pdb_chain_supfam["chain_idx"][i]):
            pdb_chain_supfam.loc[i, "chain_ipro"] = pdb_chains[pdb_chain_supfam["PDB"][i].lower()] \
                [int(pdb_chain_supfam["chain_idx"][i])]
        elif pdb_chain_supfam["chain_kegg"][i] in pdb_chains[pdb_chain_supfam["PDB"][i].lower()]:
            pdb_chain_supfam.loc[i, "chain_ipro"] = pdb_chain_supfam["chain_kegg"][i]
        else:
            print(f"kegg numbering is incorrect for {pdb_chain_supfam.PDB[i]}")
    return


def parse_supfam_local(path, pdbc=False):
    READ = False
    TAG = False
    TH = False
    rows = []
    last_field = ""
    last_tag = ""
    with open(path, "r") as f:
        for line in f:
            if "<table" in line:
                READ = True
            elif "</table>" in line:
                READ = False
            if READ:
                for c in line.strip().replace("<br>", ":"):
                    if c == "<":
                        TAG = True
                        if last_tag == "th":
                            TH = True
                        elif last_tag == "/th":
                            TH = False
                        if TH and len(last_field) > 0:
                            rows[-1].append(last_field)
                        last_tag = ""
                    elif c == ">":
                        TAG = False
                        last_field = ""
                        if last_tag == "tr":
                            rows.append([])
                    elif TAG:
                        last_tag += c
                    elif not TAG:
                        last_field += c
    if pdbc:
        pdb = []
        chain = []
    else:
        sequence = []
    supfam = []
    start = []
    end = []
    for row in rows:
        if row[0] == "Seq_ID":
            continue
        if row[3] == "-":
            print(f"row[0] is not assigned")
        for r1 in row[1].split(":"):
            if pdbc:
                r0 = row[0].split("_")
                pdb.append(r0[0].upper())
                chain.append(r0[1])
            else:
                sequence.append(row[0])
            range = r1.split("-")
            supfam.append(row[3])
            try:
                start.append(int(range[0]))
            except ValueError:
                start.append(np.nan)
            try:
                end.append(int(range[1]))
            except ValueError:
                end.append(np.nan)
    if pdbc:
        df = pd.DataFrame({"PDB": pdb, "chain_ipro": chain,
                           "Superfamily": supfam, "Start Residue": start, "End Residue": end})
    else:
        df = pd.DataFrame({"Sequence": sequence,
                           "Superfamily": supfam, "Start Residue": start, "End Residue": end})
    return df


def count_ec_supfam(df):
    sel_supfams = pd.read_csv("data/lists/supfams.csv", sep=':')
    ntp_ec = pd.read_csv("data/lists/NTP_ECs.csv")
    ec_for_supfam = {}
    for s in sel_supfams.supfam:
        ec_for_supfam[s] = []
    for i, r in df.iterrows():
        if r.EC in list(ntp_ec.ec) and r.Superfamily in list(sel_supfams.supfam):
            if r.EC not in ec_for_supfam[r.Superfamily]:
                ec_for_supfam[r.Superfamily].append(r.EC)
    for e in ec_for_supfam:
        print(e, ";".join(ec_for_supfam[e]))
    return


if __name__ == '__main__':
    ipr_per_chain = get_ipr_per_chain()
    get_ipr_corr(ipr_per_chain)
    all_supfams = pd.read_csv("data/pdb_chain_supfam.csv", index_col=0)
    sel_supfams = pd.read_csv("data/lists/supfams.csv", sep=':')
    sel_supfams_for_seq = {}
    for i, r in all_supfams.iterrows():
        if r.Superfamily in list(sel_supfams.supfam):
            if r.Sequence not in sel_supfams_for_seq:
                sel_supfams_for_seq[r.Sequence] = [r.Superfamily]
            else:
                if r.Superfamily not in sel_supfams_for_seq[r.Sequence]:
                    sel_supfams_for_seq[r.Sequence].append(r.Superfamily)
    clear_seqs = [x for x in sel_supfams_for_seq if len(sel_supfams_for_seq[x]) == 1]
    ntp_ec = pd.read_csv("data/lists/NTP_ECs.csv")
    ec_for_supfam = {}
    for s in sel_supfams.supfam:
        ec_for_supfam[s] = []
    for i, r in all_supfams.iterrows():
        if r.EC in list(ntp_ec.ec) and r.Superfamily in list(sel_supfams.supfam) and r.Sequence in clear_seqs:
            if r.EC not in ec_for_supfam[r.Superfamily]:
                ec_for_supfam[r.Superfamily].append(r.EC)
    for e in ec_for_supfam:
        print(e, ";".join(ec_for_supfam[e]))

    count_ec_supfam(all_supfams)
    ec = parse_ipr_ec()
    phosphate = pd.read_csv("data/phosphate_count.csv", index_col=0)
    repr_supfam = pd.read_csv("data/lists/repr_supfam.list").to_numpy().flatten()
    with open("list", 'r') as f:
        lines = f.readlines()
    iprs = [l.strip() for l in lines]
    for ipr in iprs:
        if not os.path.isfile(f"data/ipr/{ipr}_proc.log"):
            process_ipr(ipr, all_supfams, f"data/ipr/{ipr}_proc.log")
        parse_ipr_proc(ipr)
