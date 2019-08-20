import sys, os
from multiprocessing import Pool
import subprocess
import math
import argparse
import itertools
import pandas as pd
import scipy.stats as sc
from kpal.klib import Profile
import numpy as np
import re
from Bio import SeqIO

def get_seqs_from_fasta(sequence_path):
	records = list(SeqIO.parse(sequence_path, "fasta"))	
	seqs = []
	for x in records:
		seqs.append([x.id, str(x.seq)])
	return seqs

def get_probes_from_seq(sequence, length):
	pam = re.compile("[ATGC]GG")
	positions = [m.start() for m in pam.finditer(sequence)]
	probes = []
	for i in positions:
		if i >= length:
			probes.append(sequence[(i - length):(i + 4)])
		else:
			continue
	return probes

def calculate_GC(probe):
		GC_all = (probe.count("G") + probe.count("C")) / float(len(probe))
		GC_one = (probe[:9].count("G") + probe[:9].count("C")) / float(len(probe[:9]))
		GC_two = (probe[10:].count("G") + probe[10:].count("C")) / float(len(probe[10:]))
		return round(GC_all, 4), round(GC_one, 4), round(GC_two, 4)
		
def calc_ent(probe):
	data = pd.Series(list(probe))
	p_data= data.value_counts()/len(data)
	ent = sc.entropy(p_data)
	return ent

def kmer_complexity(probe):
	p = Profile.from_sequences([probe], 3)
	counts = p.counts
	score = np.count_nonzero(counts == 1) / (len(probe) - 2.0)
	return score

def make_eff_prediction(score, GC_all, GC_one, GC_two, ent, complexity, mismatches, coverage):
	values = [str(x) for x in [score, score, score, score, GC_all, GC_one, GC_two, complexity, ent, mismatches, coverage]]
	all = "Rscript " + "models_krispr.R " + ' '.join(values)
	#print(all)
	proc = subprocess.Popen([all], shell=True, stdout=subprocess.PIPE)
	#proc.wait()
	out = proc.stdout.read()
	out = out.decode("utf-8")
	eff = out.split('\n')[0].split()[1].rstrip().lstrip()
	if eff == 'NA':
		eff = "0"
	else:
		eff = str(round(float(eff), 4))
	
	return eff
	
def make_mutations(probe, n):
	def mutate(probe):
		p = list(probe)
		for s in range(0, 12):
			for x in ['A', 'C', 'T', 'G']:
				p[s] = x
				yield ''.join(p)
			p = list(probe)
	
	mutations = []
	if n == 3:
		for m1 in mutate(probe):
			for m2 in mutate(m1):
				for m3 in mutate(m2):
					mutations.append(m3)
	elif n == 2:
		for m1 in mutate(probe):
			for m2 in mutate(m1):
				mutations.append(m2)
	elif n == 1:
		for m1 in mutate(probe):
			mutations.append(m1)
	return list(set(mutations))

def	analyze_probes(probes, jelly_db, threads):
	tmp_file = open("tmp_probes.fasta", "w")
	
	for i in probes:
		tmp_file.write(">" + i + "\n" + i + "\n")
	tmp_file.close()
	
	all = "jellyfish " + "query -L -s tmp_probes.fasta -o tmp_out_probes.tsv " + jelly_db
	proc = subprocess.Popen([all], shell=True)
	proc.wait()
	
	out_file = open("tmp_out_probes.tsv", "r")
	results = []
	for r in out_file:
		results.append((r.split()[0], int(r.split()[1])))
	out_file.close()
	
	os.remove("tmp_probes.fasta")
	os.remove("tmp_out_probes.tsv")
	
	return results

def make_shifts(probe):
	shifts = []
	mer_len = 21
	
	nb_shifts = len(probe) - mer_len
	
	for i in range(0, nb_shifts + 1):
		shifts.append(probe[i:])
	
	return shifts

def mutate_NGG(probe):
	AGG = probe[:len(probe)-4] + "A" + probe[-3:]
	TGG = probe[:len(probe)-4] + "T" + probe[-3:]
	GGG = probe[:len(probe)-4] + "G" + probe[-3:]
	CGG = probe[:len(probe)-4] + "C" + probe[-3:]
	return AGG, TGG, GGG, CGG

def calculate_values(probe, kindex, mutations, coverage, threads):
	sum = 0
	mutated_NGG = mutate_NGG(probe)
	
	probes_to_analyze = []
	
	if mutations > 0 and mutations <= 3:
		for p in mutated_NGG:
			for m in make_mutations(p, mutations):
				probes_to_analyze.append(make_shifts(m))
				
	elif mutations == 0:
		for p in mutated_NGG:
			probes_to_analyze.append(make_shifts(p))
			
	else:
		print("mismatches not in range")
		exit(1)
	
	hits = analyze_probes(list(itertools.chain.from_iterable(probes_to_analyze)), kindex, threads)
	for h in hits:
		sum += h[1]
	
	
	corrected_sum = math.ceil(sum / coverage)
	
	GC_all, GC_one, GC_two = calculate_GC(probe)
	
	if probe[0] == "A" or probe[0] == "G":
		correct_start = "yes"
	else:
		correct_start = "no"
	
	if probe[-1] == "G":
		correct_x = "no"
	else:
		correct_x = "yes"
	
	ent = calc_ent(probe)
	
	complexity = kmer_complexity(probe)
	
	return (probe, corrected_sum, round(GC_all, 4), round(GC_one, 4), round(GC_two, 4), round(ent, 4), round(complexity, 4), correct_start, correct_x)


if __name__ == "__main__":	
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("mode", help="single: only 1 target is analyzed; multi: all targets in a fasta file are being searched and analyzed", choices=['single', 'multi'])
	parser.add_argument("-j", "--jelly_db", help="the jellyfish database file (.jf)", required = True)
	parser.add_argument("-q", "--query", help="target sequence to analyze; when mode is single: must be the target string; when mode is multi: must be a fasta file", required = True)
	parser.add_argument("-m", "--mismatches", help="number of mismatches to search", default = 3, type = int, choices=[0,1,2,3])
	parser.add_argument("-c", "--coverage", help="coverage of the jelly db", default = 1, type = int)
	parser.add_argument("--threads", help=argparse.SUPPRESS, default = 1, type = int) #"threads to use for search"
	parser.add_argument("-t", "--tabular", help="tabular output; when mode is multi, output is always tabular", action='store_true')
	parser.add_argument("-l", "--length", help="the length of the target sequence", default = 20, type = int, choices = [19,20,21])
	args = parser.parse_args()
	
	print("\nThis is KRISPR -- CRISPR offtarget and efficiency prediction with k-mers", file=sys.stderr)
	print("\nDisclaimer: These values are only predictions and can not be taken as ground truth.\n\n", file=sys.stderr)
	
	if args.mode == "single":
		if len(args.query) != args.length+4:
			print("Sequence matches not the given length")
			exit(1)
		if args.tabular == True:
			print("target\tscore\tGC_all\tGC_dist\tGC_prox\tent\tcomplexity\tstart_AG\tPAMX_G\te_eff", file=sys.stderr)
			res = [str(x) for x in calculate_values(args.query, args.jelly_db, args.mismatches, args.coverage, args.threads)]
			e_eff = make_eff_prediction(res[1], res[2], res[3], res[4], res[5], res[6], args.mismatches, args.coverage)
			print('\t'.join(res) + '\t' + str(e_eff))
		else:
			res = [str(x) for x in calculate_values(args.query, args.jelly_db, args.mismatches, args.coverage, args.threads)]
			good = ["" for _ in range(0,9)] #+
			
			eff = make_eff_prediction(res[1], res[2], res[3], res[4], res[5], res[6], args.mismatches, args.coverage)
			
			print('%-27s%-26s%-1s' % ("", "", "")) #value good
			print('%-27s%-26s%-1s' % ("target: ", res[0], ""))
			print('%-27s%-26s%-1s' % ("mismatches: ", args.mismatches, "")) #+
			print('%-27s%-26s%-1s' % ("off-target score : ", res[1], good[0]))
			print('%-27s%-26s%-1s' % ("GC all: ", res[2], good[0]))
			print('%-27s%-26s%-1s' % ("GC distal: ", res[3], good[0]))
			print('%-27s%-26s%-1s' % ("GC proximal: ", res[4], good[0]))
			print('%-27s%-26s%-1s' % ("entopy: ", res[5], good[0]))
			print('%-27s%-26s%-1s' % ("complexity: ", res[6], good[0]))
			print('%-27s%-26s%-1s' % ("starts with A or G: ", res[7], good[0]))
			print('%-27s%-26s%-1s' % ("X on PAMX ends not with G: ", res[8], good[0]))
			print("\nestimated efficiency: " + eff)
	elif args.mode == "multi":
		print("seqid\ttarget\tscore\tGC_all\tGC_dist\tGC_prox\tent\tcomplexity\tstart_AG\tPAMX_G\te_eff", file=sys.stderr)
		
		for sq in get_seqs_from_fasta(args.query):
			if len(sq[1]) < args.length+4:
				print(sq[0] + '\t' + '\t'.join(['--' for _ in range(0,7)]))
				continue
			pfs = get_probes_from_seq(sq[1], args.length)
			if not pfs:
                        	print(sq[0] + '\t' + '\t'.join(['-' for _ in range(0,7)]))
                        	continue
			for q in pfs:
				res = [str(x) for x in calculate_values(q, args.jelly_db, args.mismatches, args.coverage, args.threads)]
				e_eff = make_eff_prediction(res[1], res[2], res[3], res[4], res[5], res[6], args.mismatches, args.coverage)
				print(sq[0] + '\t' + '\t'.join(res) + '\t' + str(e_eff))
				del res
				del e_eff
	
