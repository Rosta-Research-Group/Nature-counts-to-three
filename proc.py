import argparse
import numpy as np
import mdtraj as md
from multiprocessing import Pool


def process_pdb(file):
    """
    Processes a single PDB file, can be used for trivial paralellization
    :param file: path to the PDB file
    :type file: string
    :return: void
    :rtype:
    """
    try:
        s = md.load(file)
        # select ions we like
        ions = s.top.select("resname MG or resname CA or resname MN")
        # save the pythonic residue index that is used in s.top._residues
        ionres = [s.top._atoms[x].residue.index for x in ions]
        prot = s.top.select("protein")
        found = False
        # checking every residue
        for r in s.top.residues:
            exclude = False
            p = 0
            for a in r.atoms:
                # if any of the atoms is a phosphorus, count them
                if a.element.symbol == "P":
                    p += 1
            if p > 0:
                # so this residue has some P
                pres = s.top.select("resid {0:d}".format(r.index))
                # we skip the residue names of nucleic acids (DNA/RNA)
                if s.top._atoms[pres[0]].residue.name in ["PO4", "DA", "DT", "DG", "DC", "A", "U", "G", "C"]:
                    continue
                elif p > 1:
                    # check where there is more than 1 P, select them first
                    p_atoms = s.top.select("resid {0:d} and element P".format(r.index))
                    # all this section does is checking how far P atoms are from each other to filter for
                    # phosphate anhidride chains, we don't care if they are not connected to each other
                    ppairs = []
                    for p1 in p_atoms:
                        for p2 in p_atoms:
                            ppairs.append([p1, p2])
                    pdist = md.compute_distances(s, ppairs).reshape(p, p)
                    if p > 3:
                        exclude = True
                        for l in pdist:
                            if np.sort(l)[2] <= 0.35:
                                exclude = False
                        if exclude:
                            print(f"{file} residue: {r.__str__()} does not seem to have a PPP chain, despite having {p}"
                                  f" phosphates")
                    else:
                        for l in pdist:
                            if np.sort(l)[1] > 0.35:
                                print("{1:s} residue: {0:s} 2-3 phosphates, disconnected".format(r.__str__(), file))
                                exclude = True
                    # end of connectivity check
                if not exclude:
                    # we enumerate all the ions in the structure to later check if they are close to the phosphate res
                    pairs = []
                    for i in ionres:
                        pairs.append([r.index, i])
                    # basic info of the phosphate res
                    coordination = "{2:s} p: {3:d}: Residue index: {0:d}, residue: {1:s}, ions:".format(r.index,
                                                                                                        r.__str__(),
                                                                                                        file.split(".")[0],
                                                                                                        p)
                    if len(pairs) != 0:
                        # calculates ion-phosphate res distances and adds to the print if they are close
                        contacts = md.compute_contacts(s, contacts=pairs, ignore_nonprotein=False)
                        dist = contacts[0][0]
                        for i in range(len(dist)):
                            if dist[i] < 0.5:
                                coordination += " {0:s} ({1:d}, {2:d})".format(s.top._residues[pairs[i][1]].__str__(),
                                                                               pairs[i][1],
                                                                               s.top._residues[pairs[i][1]]._atoms[0].serial)
                    # getting pythonic index for neighbouring atoms
                    nb = md.compute_neighbors(s, 0.5, pres, haystack_indices=prot)
                    cid = []
                    resid = []
                    # collecting pythonic index of chains and reisdues based on the atoms, with redundancy
                    for i in nb[0]:
                        cid.append(s.top._atoms[i].residue.chain.index)
                        # print(s.top._atoms[i].residue)
                        resid.append(s.top._atoms[i].residue.index)
                    # printing non-redundant residue names in the neigbourhood
                    for r in list(set(resid)):
                        print(s.top._residues[r])
                    # identifying the most frequent chain
                    chainid = max(set(cid), key=cid.count)
                    # getting the id of the first residue in said chain
                    chainstartid = s.top._chains[chainid]._residues[0].index
                    # this section gets the range of residue ids that are in contact of the phosphate res
                    # this was used to identify a superfamily, if the chain has multiple ones assigned to different
                    # parts
                    mincontact = 999
                    maxcontact = 0
                    for i in range(len(cid)):
                        if cid[i] == chainid:
                            if resid[i] - chainstartid < mincontact:
                                mincontact = resid[i] - chainstartid
                            if resid[i] - chainstartid > maxcontact:
                                maxcontact = resid[i] - chainstartid
                    # chain_letter = s.top._chains[chainid].chain_id # does not work in mdtraj 1.9.7
                    mincontact += 1
                    maxcontact += 1
                    # output everything
                    print("{0:s} sequence (chain {2:d} {5:.2f}): {1:s} range: {3:d} {4:d}".format(coordination,
                                                                                          s.top.to_fasta(chainid),
                                                                                          chainid,
                                                                                          mincontact, maxcontact,
                                                                                          cid.count(chainid) / len(cid)))
                    found = True
        if not found:
            print(f"{file} has no residue we care about")
    except ValueError as ve:
        print(file, "parsing problem", ve)
    except IndexError as ie:
        print(file, "parsing problem", ie)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs="+")
    args = parser.parse_args()
    pool = Pool(processes=8)
    pool.map(process_pdb, args.file)
