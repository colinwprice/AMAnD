import pandas as pd
import numpy as np
import subprocess
import time
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from taxonomy import taxid2parent

class Annoy():
    def __init__(self):

        # SET ANNOY LOCATIONS
        self.annoy_lookup = 'gen5_lookup.py'

        # SET TARGETS TO SEARCH FOR
        self.target_list = [
            {'NAME': 'AV', 'TAXID': 11019, 'SUBTAX': [59301, 374990, 37124, 371094, 59300, 97469, 11029, 11030, 11033, 177899, 59304, 11034, 84589, 78540, 11021, 1440170, 11020, 59302, 2169701, 11036, 59305, 48544, 44158, 11024, 48543, 1231903, 11023, 11039, 1159195, 332096, 1930825],'MODEL': 'alphavirus.ann', 'TAX_PICKLE': 'alphavirus_taxid_index.pickle'},
            {'NAME': 'BA', 'TAXID': 1386, 'SUBTAX': [86661, 1396, 1217984, 405535, 572264, 361100, 405532, 405534, 1003239, 1454382, 526977, 1179100, 637380, 222523, 288681, 451709, 334406, 405531, 269801, 347495, 226900, 1428, 1442, 29339, 570416, 1261129, 714359, 1279365, 180850, 29338, 1452729, 1423143, 1441, 1440, 527021, 29337, 930170, 1195464, 1217737, 1432, 1286404, 1001229, 541229, 1218175, 412694, 529122, 180856, 281309, 1392, 673518, 1449979, 743835, 1452727, 1392837, 261591, 592021, 261594, 768494, 568206, 260799, 198094, 580165, 315749, 1405, 79883, 1408, 324767, 1367477, 653685, 1938374, 1390, 692420, 1402, 279010, 1452, 1423, 135461, 224308, 1547283, 1413, 649639, 1404, 1348623, 1398, 1121088, 86664, 1471, 796606, 1664069, 632773, 1479],'MODEL': 'b_cereus_group.ann', 'TAX_PICKLE': 'b_cereus_group_taxid_index.pickle'},
            {'NAME': 'BH', 'TAXID': 32008, 'SUBTAX': [111527, 57975, 1403345, 1250339, 1249667, 1249665, 1249664, 1249663, 1249662, 1249661, 1249660, 1241582, 441157, 271848, 28450, 1306418, 1249477, 1435984, 1439852, 1435365, 557724, 1306419, 1249475, 1249658, 331978, 1249473, 1306421, 441156, 1249474, 1249476, 1439854, 1439855, 360118, 1249471, 1249472, 1249468, 1306417, 1306420, 1249659, 1335307, 320372, 536230, 320373, 884204, 357348, 1241583, 1229785, 272560, 13373, 412022, 320388, 320389, 243160, 87882, 292, 983594, 101571, 1249668, 28095],'MODEL': 'burk.ann', 'TAX_PICKLE': 'burk_taxid_index.pickle'},
            {'NAME': 'BR', 'TAXID': 234, 'SUBTAX': [29461, 204722, 1112912, 645170, 470137, 235, 359391, 1104320, 430066, 35802, 262698, 29459, 1029825, 941967, 703352, 645169, 546272, 644337, 224914],'MODEL': 'brucella.ann', 'TAX_PICKLE': 'brucella_taxid_index.pickle'},
            {'NAME': 'CB', 'TAXID': 1485, 'SUBTAX': [1509, 471871, 1075091, 378516, 1502, 195102, 195103, 1491, 1415774, 1408283, 941968, 36827, 591968, 498213, 515621, 935198, 36830, 508767, 36826, 498214, 536232, 441770, 413999, 441771, 36831, 758678, 441772, 929506, 2614128, 1581181, 1163671, 394958, 1501, 86416, 1262449, 1534, 431943, 1519, 169679, 1345695, 1493, 573061, 1561, 1415775, 1542, 386415, 36745, 1335640, 931276, 84022, 84023, 1341692, 1488, 272562, 1492, 1513, 212717, 238834, 1552, 217159, 536227, 1520, 46867, 1351755, 29341, 1497, 1216932],'MODEL': 'clostridium.ann', 'TAX_PICKLE': 'clostridium_taxid_index.pickle'},
            {'NAME': 'FT', 'TAXID': 263, 'SUBTAX': [119857, 1432652, 393011, 376619, 458234, 351581, 1232394, 264, 1450527, 676032, 1452728, 984129, 1386968, 401614, 119856, 1001534, 1001542, 543737, 177416, 1341656, 1133671, 393115, 418136, 510831, 135248, 441952], 'MODEL': 'f_tularensis.ann', 'TAX_PICKLE': 'f_tularensis_taxid_index.pickle'},
            {'NAME': 'OP', 'TAXID': 10242, 'SUBTAX': [10256, 28871, 28874, 28873, 10245, 160796, 10244, 619591, 12643, 10255, 10243],'MODEL': 'orthopoxvirus.ann','TAX_PICKLE': 'orthopox_taxid_index.pickle'},
            {'NAME': 'VC', 'TAXID': 662, 'SUBTAX': [666, 127906, 1134456, 686, 243277, 1433144, 914149, 593588, 579112, 412614, 1420885, 345073, 717610, 670, 223926, 663, 1219076, 680, 338187, 28173, 190893, 672, 45658, 687, 55601, 42323, 882102, 676, 553239, 1435069, 693153, 575788, 1891919, 29498, 1051646],'MODEL': 'vibrio.ann', 'TAX_PICKLE': 'vibrio_taxid_index.pickle'},
            {'NAME': 'YP', 'TAXID': 629, 'SUBTAX': [1649845, 632, 1345702, 1345708, 1345705, 1345706, 649716, 748678, 1345703, 1345704, 1345707, 637385, 385964, 1234662, 547048, 386656, 755858, 412420, 187410, 360102, 1035377, 637382, 1234659, 229193, 385966, 377628, 637386, 349746, 214092, 29485, 527004, 29486, 630, 150052, 393305],'MODEL': 'yersinia.ann', 'TAX_PICKLE': 'yersinia_taxid_index.pickle'},
            {'NAME': 'RT', 'TAXID': 780, 'SUBTAX': [114292, 782, 1105098, 1105094, 1290428, 449216, 1290427, 1105096, 1105099, 1105095, 1105097, 272947, 785, 257363, 114277, 783, 1105102, 1105100, 1105101, 1105103, 392021, 452659, 1105104, 1105105, 781, 272944, 42862, 315456, 33989, 787, 1105110, 786, 293614, 35790, 652620, 33992, 1105113, 1129742, 788, 293613, 33990, 336407],'MODEL': 'rickettsia.ann', 'TAX_PICKLE': 'rickettsia_taxid_index.pickle'},
            {'NAME': 'CX', 'TAXID': 776, 'SUBTAX': [777, 1293501, 360116, 360115, 434923, 434924, 434922, 227377, 2676648, 325775],'MODEL': 'coxiella.ann', 'TAX_PICKLE': 'coxiella_taxid_index.pickle'}
        ]

    def print_message(self, msg, silent, start, logfile, errorout=0):
        message = "[%s] %s\n" % (self.timeSpend(start), msg)
        with open( logfile, "a" ) as f:
            f.write( message )
            f.close()
        if errorout:
            sys.exit( message )
        elif not silent:
            sys.stderr.write( message )

    def timeSpend(self, start):
        done = time.time()
        elapsed = done - start
        return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

    # GET INPUT FILES AND SEQUENCE DATA
    # Create dataframe of ids and sequences to match with sam
    # Output only reads that match a read of interst
    #!! Phil mentioned a popen method which may be faster
    def get_input(self):

        # Build full list of taxids of interest
        full_taxids = []
        for t in self.target_list:
            full_taxids.append(t['TAXID'])
            for x in t['SUBTAX']:
                full_taxids.append(x)

        # Get Sam File
        sam_df = pd.read_csv(self.samfile,sep='\t', header=None, usecols=[0,2])
        sam_df.columns = ['read_id', 'taxid']
        sam_df['taxid'] = sam_df['taxid'].apply(self.parse_sam_column)

        # Subset Sam File on full tax list
        sam_df = sam_df[sam_df['taxid'].isin(full_taxids)]

        # Get Fastq
        ids, seqs = [], []
        for i in self.fastq:
            for seq_r in SeqIO.parse(i, 'fastq'):
                ids.append(seq_r.id)
                seqs.append(seq_r.seq)
        fastq_df = pd.DataFrame({'read_id': ids, 'seq': seqs})
        fastq_df['read_id'] = fastq_df['read_id'].apply(self.clean_reads) # DWGSim adds /1 to the read id which are not in the sam file

        # Merge Sam and Fastq
        #return pd.merge(sam_df, fastq_df, on='read_id', how='left')
        merged = pd.merge(sam_df, fastq_df, on='read_id', how='left')

        return merged.drop_duplicates()

    # PULL ALL CHILD TAXIDS  --  need to implement into tax-library
    # Using rollup.py and sqlite database, grab all sub genus taxids of a tree
    #def get_taxids(self, taxid):
    #    query = FromTaxonomy()
    #    try:
    #        query.get(int(taxid))
    #        return query.childlist(int(taxid))
    #    except:
    #        print('{} not in taxonomy database.'.format(taxid))

    def parse_sam_column(self, sam_input):
        """
        # This is just a quick function to clean up sam columns
        :param input:
        :return taxid:
        """
        sam_input = sam_input.split('|')
        if len(sam_input) > 1:
            output = int(float(sam_input[-2]))
        else:
            output = input[0]
        return output

    def clean_reads(self, read_id):
        return read_id.split('/')[0]

    # GET PARENT TAXID
    # If taxid not in report tsv lookup closest parent tax that is in reports
    def get_parent(self, this_tax, results_tax_list):
        # Remove the float if it has one
        this_tax = int(this_tax)

        # See if the strain tax decimal <taxid.1> is messing up matching
        if this_tax in results_tax_list:
            return this_tax

        # get parent of current tax, if not in results get keep looping
        parent_tax = False
        while parent_tax == False:
            try:
                next_tax = taxid2parent(this_tax)
            except:
                parent_tax = True
                return this_tax
            #print('Getting Parent Taxid: \t{}\t\t=\t{}'.format(this_tax, next_tax))  -  I was getting stuck in a loop
            if next_tax in results_tax_list:
                parent_tax = True
                return next_tax
            # To prevent non-stop loops if there is an error
            if next_tax == 0 or this_tax == 0 or next_tax == this_tax:
                parent_tax = True
                return next_tax
            this_tax = next_tax

    # START THE ANNOY MODEL
    # This is the helper main function of annoy that outlines the steps of the annoy pipeline
    # Will contain the variable and also output the results
    def Start(self):

        self.print_message( "Beginning ANNOY analysis of PanGIA results.", self.silent, self.begin_t, self.logfile )

        # Parse Input Files
        # Returns df with [read_id, taxid, seq]
        self.print_message( 'Parsing fastq(s) and sam file.', self.silent, self.begin_t, self.logfile )
        sam_df = self.get_input()

        annoy_results = pd.DataFrame()

        total_num_annoy_reads = 0
        total_num_deleted_reads = 0

        # Loop through each target
        for target in self.target_list:
            # Get child taxids of target
            tax_list = target['SUBTAX']
            tax_list.append(target['TAXID'])

            # Subset dataframe of reads for taxids associated to each target
            this_reads = sam_df[sam_df['taxid'].isin(tax_list)]

            # If there are no reads for this target go to the next target
            if len(this_reads) == 0:
                self.print_message( 'No reads found for {}.'.format(target['NAME']), self.silent, self.begin_t, self.logfile )
                continue

            self.print_message( 'Found reads specific to {}.'.format(target['NAME']), self.silent, self.begin_t, self.logfile )

            # Remove reads that are too short
            num_of_treads = len(this_reads)
            this_reads = this_reads[this_reads['seq'].apply(lambda x: len(x)>123)]

            # Calculate how many reads are deleted. pass that info to user.
            # If no reads left don't run annoy.
            deleted_reads = num_of_treads - len(this_reads)
            total_num_annoy_reads += num_of_treads
            total_num_deleted_reads += deleted_reads

            if len(this_reads) == 0:
                self.print_message( '\t All {} reads removed due to length.'.format(num_of_treads), self.silent, self.begin_t, self.logfile )
                continue
            else:
                if deleted_reads == 0:
                    text_del_percent = 0.0
                else:
                    text_del_percent = round(int(num_of_treads)/int(deleted_reads)*100, 1)
                self.print_message( '\t{}({}%) {} reads of {} removed due to length.'.format(
                    deleted_reads,
                    text_del_percent,
                    target['NAME'],
                    num_of_treads
                ), self.silent, self.begin_t, self.logfile )

            # Create fasta file in tmp directory
            this_fasta = []
            for index, row in this_reads.iterrows():
                this_fasta.append(SeqRecord(Seq(str(row['seq'])), id=str(row['read_id']), name=str(row['read_id']), description='', dbxrefs=[]))
            tmp_fasta = '{}/{}.fasta'.format(self.tmpdir, target['NAME']) # fastq tmp file name
            SeqIO.write(this_fasta, tmp_fasta, 'fasta')

            self.print_message( '\tRunning annoy index for {} on {} reads.'.format(target['NAME'], len(this_reads)), self.silent, self.begin_t, self.logfile )

            # Run Annoy
            tmp_annoy_tsv = '{}/{}_annoy.tsv'.format(self.tmpdir, target['NAME'], '') # annoy results file name
            annoy_cmd = 'python {} -a {} -t {} -r {} -o {}'.format(
                '{}/{}'.format(Path(__file__).parent.absolute(), self.annoy_lookup),
                '{}annoy/{}'.format(self.dbPath, target['MODEL']),
                '{}annoy/{}'.format(self.dbPath, target['TAX_PICKLE']),
                tmp_fasta,
                tmp_annoy_tsv
            )

            # Annoy subprocess
            # Outputs temp file (annoy_tsv)
            annoy_subprocess = subprocess.check_output([annoy_cmd], shell=True)
            annoy_df = pd.read_csv(tmp_annoy_tsv, sep='\t')
            this_counter = Counter(annoy_df['taxid'].tolist())
            counter_tax, counter_num = [], []
            for x, y in this_counter.items():
                counter_tax.append(x)
                counter_num.append(y)

            this_results = pd.DataFrame({'taxonomy':counter_tax,'read_count':counter_num})
            # Concat results
            annoy_results = pd.concat([annoy_results,this_results])

        # Sort annoy results and combine read_count for any duplicate taxids
        results_df = pd.read_csv(self.outfile, sep='\t')
        results_df['ANNOY_RC'] = int(0)
        results_df['ANNOY'] = ""

        # Only continue if there are results for Annoy
        if len(annoy_results) > 0:

            annoy_sum = annoy_results.groupby(['taxonomy'])['read_count'].sum().to_frame()

            # Output annoy classification tsv
            annoy_tsv_file = str(self.outfile).replace('.report.tsv', '.classification.tsv')
            annoy_tsv = annoy_sum.reset_index()
            annoy_tsv.to_csv(annoy_tsv_file, sep='\t', index=False)

            # Display info about the run to screen
            if total_num_deleted_reads == 0:
                text_del_percent = 0.0
            else:
                text_del_percent = round(int(total_num_annoy_reads)/int(total_num_deleted_reads)*100, 1)
            text = 'Adding annoy results to output report file.\nNumber of reads for target taxids:\t {}\nNumber of reads too short to process:\t {} ({}%)\n\n-- ANNOY RESULTS --\n{}\n'.format(
                total_num_annoy_reads,
                total_num_deleted_reads,
                text_del_percent,
                annoy_tsv
            )
            self.print_message( text, self.silent, self.begin_t, self.logfile )

            # Combine annoy tsv and results tsv
            annoy_sum.columns = ['ANNOY_RC']
            annoy_sum = annoy_sum.reset_index()
            # Turn taxonomy to numeric if string turn to nan
            annoy_sum['taxonomy'] = pd.to_numeric(annoy_sum['taxonomy'], errors='coerce')
            # Remove rows with taxonomy = nan
            annoy_sum = annoy_sum[annoy_sum['taxonomy'].notna()]
            results_df['taxonomy'] = results_df['TAXID']

            for i, r in annoy_sum.iterrows():
                if r['taxonomy'] == np.nan:
                    continue
                this_taxid = r['taxonomy']

                # Remove .0 decimals
                if this_taxid % 1 == 0:
                    this_taxid = int(this_taxid)

                if len(results_df[results_df['taxonomy'] == this_taxid]) > 0:
                    results_df.loc[results_df['taxonomy'] == this_taxid, 'ANNOY_RC'] += int(r['ANNOY_RC'])
                    results_df.loc[results_df['taxonomy'] == this_taxid, 'ANNOY'] += '{}:{};'.format(this_taxid, int(r['ANNOY_RC']))
                else:
                    parent_tax = self.get_parent(r['taxonomy'], results_df['taxonomy'].tolist())
                    if parent_tax != 0:
                        results_df.loc[results_df['taxonomy'] == parent_tax, 'ANNOY_RC'] += int(r['ANNOY_RC'])
                        results_df.loc[results_df['taxonomy'] == parent_tax, 'ANNOY'] += '{}:{};'.format(this_taxid, int(r['ANNOY_RC']))
                    else:
                        print('Failed to find parent for taxid:{}'.format(this_taxid))

            del results_df['taxonomy']

        # There were no ANNOY results
        else:
            if total_num_deleted_reads == 0:
                text_del_percent = 0.0
            else:
                text_del_percent = round(int(total_num_annoy_reads)/int(total_num_deleted_reads)*100, 1)

            text = 'No ANNOY results. {} of {} ({}%) reads removed due to length.'.format(
                total_num_deleted_reads,
                total_num_annoy_reads,
                text_del_percent
            )
            self.print_message(text, self.silent, self.begin_t, self.logfile )

        # Write out new report.tsv
        results_df.to_csv(self.outfile, sep='\t', index=False)
