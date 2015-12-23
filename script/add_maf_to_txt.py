#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: add-maf-to-txt.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Mon 13 Apr 2015 04:48:21 PM CST
################################################


import os,sys
'''input anno txt file'''

##dbsnp reference
dbsnp = '/project/liucj/REFDATA/dbSNP/hg19/snp142.txt'
##index of dbsnp
index = '/project/liucj/REFDATA/dbSNP/hg19/snp142.txt.idx.sort'


def change(ob):
	tmp = list()
	alpha = {
			'A':'T',
			'G':'C',
			'C':'G',
			'T':'A'
			}
	for i in ob:
		if i in alpha:
			tmp.append(alpha[i])
		else:
			tmp.append(i)
	
	tmp = ''.join(tmp)
	return tmp
		
def convert(rec):
	'''change the pos and strand'''
	if rec[11] != 'insertion':
		rec[2] = str(int(rec[2]) + 1)

	if rec[6] == '-':
		rec[9] = change(rec[9])
		if int(rec[21]) > 0:
			rec[22] = change(rec[22])
	
	return rec

def get_freq(freq):
	al = freq[0].split(',')
	ct = freq[1].split(',')
	fq = freq[2].split(',')

	result = dict()

	for i in range(len(al)):
		result[al[i]] = [ct[i], fq[i]]
	
	return result

def proc_rec(arr, rec):
	'''use refUCSC'''
	result = ['.' for i in range(5)]

	flag = 0

	rec = convert(rec)
	ob = rec[9].split('/')

	if arr[0:3] == rec[1:4] and arr[3] == rec[8] and arr[4] in ob:
		flag = 1
		result[0] = rec[4]
		if rec[21] != '0':
			#print len(rec)
			freq = get_freq(rec[22:25])
			if arr[3] in freq:
				result[1:3] = freq[arr[3]]
			if arr[4] in freq:
				result[3:5] = freq[arr[4]]

	return flag, result

def run(f, idx = index, dbsnp = dbsnp):
	#print "Be sure of chromosome sort from 1-22,chrx,chry,chrm"
	#open file	
	f_h = open(f, 'r')
	index = open(idx, 'r')
	dbsnp = open(dbsnp, 'r')

	#initiate the index_line
	index_line = index.readline().rstrip()
	index_arr = index_line.split("\t")
	
	for line in f_h:
		if line.startswith('Chr'): continue
		
		line = line.rstrip()
		arr = line.split("\t")

		my_freq = list()
		
		#read through the idx file
		while index_arr[0] != arr[0] or (int(arr[1]) - int(index_arr[1])) > 1000:
			#print index_arr[0]
			index_line = index.readline().rstrip()
			index_arr = index_line.split("\t")

		#print index_line
		seek1 = int(index_arr[2])
		seek2 = int(index_arr[3])

		dbsnp.seek(seek1)

		while True:
			if dbsnp.tell() > seek2: break

			rec = dbsnp.readline()
			rec = rec.rstrip()
			rec_arr = rec.split("\t")
			
			flag, my_freq = proc_rec(arr, rec_arr)
			if flag == 1: break
		
		#print arr[0:5],my_freq
		
		arr.append("\t".join(my_freq))
		line = "\t".join(arr)
		print line


	#close handle
	dbsnp.close()
	index.close()
	f_h.close()

if __name__ == '__main__':
	run(sys.argv[1]) #The argv[1] is the multianno file.
	


