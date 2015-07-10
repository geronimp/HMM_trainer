import argparse, os, sys, subprocess

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
        return mapped_reads
    
    def _srch(self, hmm, seqs):
        
        cmd = 'orfm -m 96 %s | hmmsearch -o /dev/null --domtblout /dev/stdout \
               %s %s' % (hmm, seqs)
        
        hmmsearch_output = subprocess.check_output(
                                                   cmd,
                                                   shell=True
                                                   )
        
        found_read_names = [x.split('\t')[0] for x in samtools_output.split('\n') \
                            if not x.startswith('#')]
        
        import IPython ; IPython.embed()
        return found_read_names
        
    def train(self, hmm, genome_ids):
        
        for gene_id, genome_id in self._get(genome_ids):
            
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
            
            known_reads = self._freadsin(gene_id, genome_file, gff_file, bam_file)
            found_reads = self._srch(hmm, mock_file)
            
        
HMMTrainer().train(sys.argv[1], sys.argv[2])