#!/usr/bin/python

import sys
import gzip
#import xml.etree.ElementTree as ET
from xml.etree import cElementTree

odat_path = 'output/uniprot_data_mouse.tsv'
with open(odat_path, 'w') as odat:
	odat.write('acc\tgene\tlocation\ttopology\tboth\n')

in_path = sys.argv[1] # you need the file uniprot_sprot_rodents.xml.gz 

namespace = '{http://uniprot.org/uniprot}'

gene = ''

e = sys.stderr

entries_read = 0
rows_written = 0

f = gzip.open(in_path, mode='rb') # 'rb' ?
iterparser = iter(cElementTree.iterparse(f, events=("start", "end")))
for event, elem in iterparser:
	if elem.tag == namespace + 'uniprot':
		root = elem
		continue
	if elem.tag == namespace + 'entry' and event == 'end': # wait until that node is done parsing
		gene = ''
		comments = ''
		entry = elem
		entries_read += 1
		accession_number = entry.find(namespace + 'accession').text
		# only parse the human protein entries
		is_human_entry = False
		organism_node = entry.find(namespace + 'organism')
		for organism_name_node in organism_node.findall(namespace + 'name'):
			if organism_name_node.text == 'Mouse':
				is_human_entry = True
		if not is_human_entry:
			continue
		# grab gene symbol
		gene_node = entry.find(namespace + 'gene')
		if gene_node is None:
			continue
		for gene_name_node in gene_node.findall(namespace + 'name'):
			if gene_name_node.attrib['type'] == 'primary':
				gene = gene_name_node.text
		# grab location data
		location_text = []
		topology_text = []
		for comment_node in entry.findall(namespace + 'comment'):
			if comment_node.attrib['type'] == 'subcellular location':
				subloc_nodes = comment_node.findall(namespace + 'subcellularLocation')
				for subloc_node in subloc_nodes:
					location_text_nodes = subloc_node.findall(namespace + 'location')
					for location_text_node in location_text_nodes:
						location_text.append(location_text_node.text)
					topology_text_nodes = subloc_node.findall(namespace + 'topology')
					for topology_text_node in topology_text_nodes:
						topology_text.append(topology_text_node.text)
		location_string = ','.join(set(location_text))
		topology_string = ','.join(set(topology_text))
		both_string = ','.join(set(location_text + topology_text))
		with open(odat_path, 'a') as odat:
			odat.write('\t'.join([accession_number, gene, location_string, topology_string, both_string]) + '\n')
		rows_written += 1
		e.write('\r' + str(entries_read) + ' entries read, ' + str(rows_written) + ' rows written, currently on ' + gene + '...')
		elem.clear()
		root.clear()



