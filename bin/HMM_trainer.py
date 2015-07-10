import argparse, os, sys, subprocess, regex, re

PATH_TO_GENOMES = "/srv/db/img/latest/genomes/"
PATH_TO_BAMS    = "/srv/projects/abisko/annotation/joel/01_training_set/bams/"
PATH_TO_MOCKS   = "/srv/projects/abisko/annotation/joel/01_training_set/mocks/"

class HMMTrainer:
    
    def __init__(self):pass
    
    def _get(self, path):
        for id in [x.strip().split() for x in open(path, 'r').readlines()]:
            yield id[0], id[1]
    
    def _rgff(self, gff_file):
        genes_id_to_region={}
        for line in open(gff_file).readlines():
            if not line.startswith('#'):
                line=line.split('\t')
                c = line[0]
                f = line[3]
                t = line[4]
                i = line[-1].split(';')[0][3:]
                genes_id_to_region[i]=[c,f,t]
        return genes_id_to_region
        
    def _freadsin(self, gene_id, genome_file, gff_file, bam_file):
        genome_dict   = self._rgff(gff_file)
        gene_location = genome_dict[gene_id]
        cmd = '''samtools view %s "%s:%s-%s" ''' % (
                                                    bam_file, 
                                                    gene_location[0], 
                                                    gene_location[1], 
                                                    gene_location[2]
                                                    )
        
        samtools_output = subprocess.check_output(
                                                  cmd,
                                                  shell=True
                                                  )
        
        mapped_read_names = [x.split('\t')[0] for x in samtools_output.split('\n') if x]
        return mapped_read_names
    
    def _srch(self, hmm, seqs):
        found_read_names_hash = {}
        orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
        cmd = 'orfm -m 96 %s | hmmsearch -o /dev/null --domtblout /dev/stdout \
               %s -' % (seqs, hmm)
        
        
        hmmsearch_output = subprocess.check_output(
                                                   cmd,
                                                   shell=True
                                                   )
        
        for line in hmmsearch_output.split('\n'):
            if line:
                if not line.startswith('#'):
                    name = orfm_regex.match(line.split()[0]).groups(0)[0] 
                    bit  = line.split()[7]
                    if name in found_read_names_hash:
                        found_read_names_hash[name].append(float(bit))
                    else:
                        found_read_names_hash[name] = [float(bit)]
        
        found_read_names_list = [[key, max(item)] for key, item \
                                 in found_read_names_hash.iteritems()]

        return found_read_names_list
        
    def train(self, hmm, genome_ids):
        
        max_known_hits      = []
        max_false_positives = []
        false_negatives     = []
        
        for gene_id, genome_id in self._get(genome_ids):
            known_hits      = []
            false_positives = []
            genome_file = os.path.join(PATH_TO_GENOMES, 
                                       genome_id, 
                                       genome_id + '.fna')
            gff_file    = os.path.join(PATH_TO_GENOMES, 
                                       genome_id, 
                                       genome_id + '.gff')
            bam_file    = os.path.join(PATH_TO_BAMS, 
                                       genome_id + '.bam')
            mock_file   = os.path.join(PATH_TO_MOCKS,
                                       genome_id + '_sammy_mock.fna')
            
            known_reads = set(
                              self._freadsin(gene_id, 
                                             genome_file, 
                                             gff_file, 
                                             bam_file)
                              )
            
            found_reads = self._srch(
                                     hmm, 
                                     mock_file
                                     )

            kr = [x for x in known_reads]
            
            for hit in found_reads:
                hit_name = hit[0]
                hit_bit  = hit[1]
                
                if hit_name in known_reads:
                    known_hits.append(hit_bit)
                    kr.remove(hit_name)
                else:
                    false_positives.append(hit_bit)

            false_negatives = kr
            
            if known_hits:
                max_known_hits.append(min(known_hits))
            if false_positives:
                max_false_positives.append(max(false_positives))
            
        import IPython ; IPython.embed()
            
            
        
HMMTrainer().train(sys.argv[1], sys.argv[2])