#Little test for git 
import json
from Bio.Restriction import * 
from dnachisel import *
from Bio import Seq


###########################################################
# To recode refactored zone (and corresponding regions)
###########################################################

def recode_max_editing_distance(seq_ref,seq_to_change,codon_table_mesoL1,seed):
    """
    Some duplications may occur during the refactoring process of M. florum. This function help solve these problems.
    the length(sequence)%3 must be 0, otherwise the sequence won't be recoded. Sequence in frame ONLY!
    The input sequences must be the same length.
    
    This function do these steps in this priority order:
    1) Choose the most different codon which encodes for the same amino acid (larger hamming distance);
    2) Choose the codon based on GC content. The highest GC content codon is selected;
    3) Choose based on the codon frequency in mesoplasma florum;
    4) Choose randomly in remaining codons.
    
    args:
    codon_table_mesoL1 is a dictionnary in this format: {'Q': {'CAA': 0.95, 'CAG': 0.05}, 'N': {'AAT': 0.77, 'AAC': 0.23}, 'K': {'AAA'...}...}
    seq_ref is a string of DNA. It has to diverge from the output sequence has described above.
    seq_to_change is a string of DNA. It has to keep the same amino acid than the output sequence
    
    return:
    a recoded DNA string
    
    """
    #important libraries for this function
    
    from Bio.SeqUtils import GC

    import numpy as np
    
#     Set the seed :

    seq_ref_uC = seq_ref.upper()
    seq_to_change_uC = seq_to_change.upper()

    
    #important conditions for users:
    assert (len(seq_ref_uC)%3 == 0),"The length of the input sequence must be a multiple of three in an ORF"
    
    #Are sequences the same length?:
    assert (len(seq_ref_uC) == len(seq_to_change_uC)),"The input sequences must be the same length"
    
    
    #amino acid:[list of associated codons]
    bon_dico = {cle:list(valeur.keys()) for cle, valeur in codon_table_mesoL1.items()}

    # codon : associated amino acid
    reverse_dico = {codons:cle for cle, valeur in bon_dico.items() for codons in valeur}


    #codon:occurence in CDS for mesoplasma florum
    codon_occurence = {codons:occurence for acide_amine, codons_occurence in codon_table_mesoL1.items() for codons, occurence in codons_occurence.items()}

    hamming_sequence_changed = ''
    
    #scanning in the input string

    for idx in range(int(len(seq_to_change_uC)/3)):
        if seq_to_change_uC[idx*3:idx*3+3] != 'TAG':
            hamming = {}

            #for all the genes which encodes for the same amino acid:
            for codons in bon_dico.get(reverse_dico.get(str(seq_to_change_uC[idx*3:idx*3+3]))):
                score = 0
                #count the difference between each nucleotide
                for compteur in range(len(seq_ref_uC[idx*3:idx*3+3])):
                    if codons[compteur] != seq_ref_uC[idx*3:idx*3+3][compteur]:
                        score = score + 1
                    hamming.update({str(codons) : score})
            #print(hamming)
            
                #Conditions
                #Find the most different codon which encodes for the same amino acid (larger hamming distance):
            liste_distance = []
            [liste_distance.append(triplet) for triplet, differences in hamming.items() if differences == max(hamming.values())]


            if len(liste_distance) == 1 :
                choix = liste_distance[0]

                # Choose the codon based on GC content: the more it is the better it is :
            elif len(liste_distance) > 1 :
                Gc_content = {triplet : GC(triplet) for triplet in liste_distance}
                liste_GC = []
                [liste_GC.append(triplet) for triplet, gc in Gc_content.items() if gc == max(Gc_content.values())]

                if len(liste_GC) == 1:
                    choix = liste_GC[0]

                #Choose based on the codon frequency in mesoplasma florum:    
                elif len(liste_GC) > 1:
                    freq_meso = {triplet : codon_occurence.get(str(triplet)) for triplet in liste_GC}
                    liste_occu = []
                    [liste_occu.append(triplet) for triplet, occurence in freq_meso.items() if occurence == max(freq_meso.values())]

                    if len(liste_occu) == 1:
                        choix = liste_occu[0]

                    #multiple codon candidate still left (really rare case, but possible): choose randomly in the list:
                    elif len(liste_occu) > 1:
                        np.random.seed(seed+idx)
                        choix = np.random.choice(liste_occu)

        elif seq_to_change_uC[idx*3:idx*3+3] == 'TAG':
            choix = 'TAG'
            
        #We dont want a CGG in the recoded sequence of M. florum:
        if choix == 'CGG':
            choix = 'AGA'
            
        #We dont wnat TGA either (really rare case, but still):
        if choix == "TGA":
            choix = "TGG"
            
        #Add choice in the hamming sequence
        hamming_sequence_changed = hamming_sequence_changed + choix


    return hamming_sequence_changed


def recode_max_edit_left_refactored_zone(genome,recoded_genome,codon_table_meso,left_end,left_gene,right_gene,right_gene_recoded_CDS,genes_overlapping,seed):
    # here the only case involved are ++ and +-, so left gene is on + strand
    reed_zone_length = left_gene['end'] - left_end
    # adjust zone for it to fit with codons (modulo3)
    TAG_part = ''
    if reed_zone_length%3 != 0:
        length_to_mod3 = 3-(reed_zone_length%3)
        TAG_part = 'TAG'[-length_to_mod3:]
        reed_zone_length += length_to_mod3
    if genes_overlapping:
        overlap_length = left_gene['end'] - right_gene['start']
        wt_genome_part = genome[left_end:right_gene['start']]
        recoded_right_gene_part = right_gene_recoded_CDS[:overlap_length]
        ref_sequence = TAG_part + wt_genome_part + recoded_right_gene_part
    else:
        ref_sequence = TAG_part + genome[left_end:left_gene['end']]
    sequence_to_recode = recoded_genome[-(reed_zone_length):]
    recoded_sequence = recode_max_editing_distance(ref_sequence, sequence_to_recode, codon_table_meso,seed)
    recoded_genome = recoded_genome[:-(reed_zone_length)] + recoded_sequence
        
    return reed_zone_length, recoded_genome

def recode_max_edit_right_refactored_zone(genome,left_gene,left_gene_recoded_CDS,right_gene_recoded_CDS,codon_table_meso,right_end,right_gene,genes_overlapping,seed):
    # here the only case involved are -- and +-, so right gene is on - strand
    reed_zone_length = right_end - right_gene['start']
    # adjust zone for it to fit with codons (modulo3,seed
    TAG_part = ''
    if reed_zone_length%3 != 0:
        length_to_mod3 = 3-(reed_zone_length%3)
        TAG_part = Seq.Seq('TAG'[-length_to_mod3:]).reverse_complement()
        reed_zone_length += length_to_mod3
    if genes_overlapping:
        overlap_length = left_gene['end'] - right_gene['start']
        wt_genome_part = genome[left_gene['end']:right_end]
        recoded_left_gene_part = left_gene_recoded_CDS[-overlap_length:]
        ref_sequence = recoded_left_gene_part + wt_genome_part + TAG_part
    else:
        ref_sequence = genome[right_gene['start']:right_end] + TAG_part
    # need to reverse complement the ref_sequence because the right_gene is on - strand
    ref_sequence = ref_sequence.reverse_complement()
    sequence_to_recode = right_gene_recoded_CDS[-reed_zone_length:]
    recoded_sequence = recode_max_editing_distance(ref_sequence, sequence_to_recode, codon_table_meso,seed)
    return reed_zone_length, recoded_sequence


###########################################################
# To find intergenic region limits/ends
###########################################################

def find_overlapping_promoters_start(TU_df,TUs,gene,strand):
    start_promoters = []
    if strand == '+':
        for TU in TUs:
            TU = TU_df.loc[TU_df['TU name'] == TU]
            if int(TU['TU start']) -75 > int(gene['start']) and int(TU['TU start']) -75 < int(gene['end']):
                start_promoters.append(int(TU['TU start']) -75)
    else:
        for TU in TUs:
            TU = TU_df.loc[TU_df['TU name'] == TU]
            if int(TU['TU end']) +75 > int(gene['start']) and int(TU['TU end']) +75 < int(gene['end']):
                start_promoters.append(int(TU['TU end']) +75)
    return start_promoters

def get_end_promoter(TU_df,gene_to_TU,gene_possibly_overlapped,gene_with_possible_promoter):
    if gene_with_possible_promoter['strand'] == '+':
        end_promoters = [gene_possibly_overlapped['end']]
        gene_have_promoters = gene_with_possible_promoter['gene name (RAST)'] in gene_to_TU.keys()
        if gene_have_promoters:
            TUs = gene_to_TU[gene_with_possible_promoter['gene name (RAST)']]
            end_promoters = end_promoters + find_overlapping_promoters_start(TU_df,TUs,gene_possibly_overlapped,'+')
        end_promoter = min(end_promoters)
    else:
        end_promoters = [gene_possibly_overlapped['start']]
        gene_have_promoters = gene_with_possible_promoter['gene name (RAST)'] in gene_to_TU.keys()
        if gene_have_promoters:
            TUs = gene_to_TU[gene_with_possible_promoter['gene name (RAST)']]
            end_promoters = end_promoters + find_overlapping_promoters_start(TU_df,TUs,gene_possibly_overlapped,'-')
        if gene_with_possible_promoter['gene name (RAST)'] == 'peg.55':
            end_promoter = end_promoters[1]
        else:
            end_promoter = max(end_promoters)
    return end_promoter

def find_self_overlapping_terminator_start(terminator_df,gene):
    terminator = terminator_df.loc[(terminator_df['upstream gene']==gene['gene name (RAST)'])]
    if gene['strand'] == '+':
        start_terminator = min(int(gene['end']),int(terminator['term start']))
    elif gene['strand'] == '-':
        start_terminator = max(int(gene['start']),int(terminator['term stop']))
    return start_terminator

def get_end_terminator(terminator_df,gene):
    gene_have_terminator = gene['gene name (RAST)'] in list(terminator_df['upstream gene'])
    if gene_have_terminator:
        end_terminator = find_self_overlapping_terminator_start(terminator_df,gene)
    else:
        if gene['strand'] == '+':
            end_terminator = gene['end']
        else:
            end_terminator = gene['start']
    return end_terminator

###########################################################
# For Annotation files
###########################################################


def df2gff3(df, output_file_dir):
    
    """
        This function transforms a dataframe into gff3 file.
        WARNING: None of these input caracteristics are verified: garbage in = garbage out. 

        ARGS:

        a) Data frame of this format (pandas data frame);
        start     end strand gene name (RAST)
    0         0    1332      +            peg.1
    1      1559    2681      +            peg.2
    2      2707    3682      -            peg.3
    3      3750    4278      +            peg.4
    4      4279    5083      +            peg.5
    ..      ...     ...    ...              ...
    715  789718  790147      -          peg.681
    716  790280  791081      -          peg.682
    717  791180  792428      -          peg.683
    718  792444  792777      -          peg.684
    719  792789  792924      -          peg.685

    [720 rows x 4 columns]

        b) a path where the output is written (string). 


        RETURN:
        None. An output gff3 file is written during the process. 


    """

    gff3 = df.copy()
    gff3["start"] = gff3["start"]+1
    gff3["annotation"] = 'ID=gene-'+gff3['gene name (RAST)']+';Name='+gff3['gene name (RAST)']+';gbkey=Gene;gene_biotype=protein_coding;locus_tag='+gff3['gene name (RAST)']
    gff3["chromosome"] = "L1"
    gff3["annot_type"] = "RAST"
    gff3["gene"] = "gene"
    gff3["point"] = "."

    gff3 = gff3[['chromosome', 'annot_type', 'gene', 'start', 'end', 'point', 'strand', 'point', 'annotation']]

    with open(str(output_file_dir), 'w') as f_in:
        f_in.write('##gff-version 3'+'\n')
        f_in.write(gff3.to_csv(header = False, index = False, sep = '\t'))


def annotate_cds(annotations_df,gene,position):
    annotations_df = annotations_df.append({
        "start":position,
        "end":position+(gene["end"]-gene["start"]),
        "gene name (RAST)":gene["gene name (RAST)"]+"_rec",
        "strand":gene["strand"]}, ignore_index=True)
    return annotations_df

def annotate_IR(annotations_df,left_gene,left_end,right_gene,right_end,position):
    if left_gene['end'] < right_gene['start']:
        if left_end < left_gene['end']:
            annotations_df = annotations_df.append({
                "start":position,
                "end":position+(left_gene['end']-left_end),
                "gene name (RAST)":left_gene["gene name (RAST)"]+"_e_refactored_IR",
                "strand":left_gene["strand"]}, ignore_index=True)
        if right_end > right_gene['start']:
            annotations_df = annotations_df.append({
                "start":position+(right_end-left_end)-(right_end-right_gene['start']),
                "end":position+(right_end-left_end),
                "gene name (RAST)":right_gene["gene name (RAST)"]+"_s_refactored_IR",
                "strand":right_gene["strand"]}, ignore_index=True)
    else:
        if left_end < left_gene['end']:
            annotations_df = annotations_df.append({
                "start":position,
                #Warning+20?
                "end":position+20,
                "gene name (RAST)":left_gene["gene name (RAST)"]+"_e_refactored_IR",
                "strand":left_gene["strand"]}, ignore_index=True)
        if right_end > right_gene['start']:
            #"start":position+(right_end-left_end)-(right_end-right_gene['start'])
            annotations_df = annotations_df.append({
                "start":position,
                "end":position+(right_end-left_end),
                "gene name (RAST)":right_gene["gene name (RAST)"]+"_s_refactored_IR",
                "strand":right_gene["strand"]}, ignore_index=True)
    return annotations_df

def annotate_IR_wtz(annotations_df,left_gene,left_end,right_gene,right_end,position):
    if left_end < left_gene['end']:
        annotations_df = annotations_df.append({
            "start":position-(left_gene['end']-left_end),
            "end":position,
            "gene name (RAST)":left_gene["gene name (RAST)"]+"_e_wtz_IR",
            "strand":left_gene["strand"]}, ignore_index=True)

    if right_end > right_gene['start']:
        annotations_df = annotations_df.append({
            "start":position+(right_gene['start']-left_gene['end']),
            "end":position+(right_gene['start']-left_gene['end'])+(right_end-right_gene['start']),
            "gene name (RAST)":right_gene["gene name (RAST)"]+"_s_wtz_IR",
            "strand":right_gene["strand"]}, ignore_index=True)
    return annotations_df

def annotate_reedited_sequence(annotations_df,left_gene,left_end,right_gene,right_end,zone_length,position,side):
    if side == 'left':
        annotations_df = annotations_df.append({
            "start":position-zone_length,
            "end":position,
            "gene name (RAST)":left_gene["gene name (RAST)"]+"_reed",
            "strand":left_gene["strand"]}, ignore_index=True)
    else:
        annotations_df = annotations_df.append({
            "start":position+(right_end-left_end),
            "end":position+(right_end-left_end)+zone_length,
            "gene name (RAST)":right_gene["gene name (RAST)"]+"_reed",
            "strand":right_gene["strand"]}, ignore_index=True)
                                 
    return annotations_df

#Some tiny functions:
def create_rev_codon_table(all_CDSs, number_table):
    return reverse_codon_table(create_codon_table(all_CDSs, number_table))


def create_file_codon_table(codon_table,name):
    with open('../data/'+str(name), 'w') as f:
        json.dump(codon_table, f)
        
        
#Recode functions:
def codonoptimize_match_codon_usage_no_TGA_only_ATG(seq,initial_codon_usage,target_codon_usage,seed):
    
    """
    This function optimize recode the M.florum CDSs based on the codon table of B.subtilis. It removes homopolymers and other features present in WT CDSs sequences.  
    
    ARGS:
    One CDS at a time (string)
    Initial_codon_usage (in our case, M.florum reverse codon table)
    Target_codon_usage (in our cas b. Subtilis reverse codon table)
    idx = for the seed and randomness of DNA chisel. 
    
    RETURN:
    Recoded CDS sequence of M.florum
    
    """
    import sys
    import numpy as np
    from Bio import Seq
    
#     Set the seed :
    np.random.seed(seed)
  
    #To remove the sdterr (really annoying)
    orig_stderr = sys.stderr
    f = open('out.txt', 'w')
    sys.stderr = f
    #To remove the sdterr (really annoying)
    
    no_TGA_seq = ""
    for i in range(int(len(seq)/3)):
        if i == 0:
            no_TGA_seq = "ATG"
            continue
        if seq[(i*3):(i*3)+3] == "TGA":
            no_TGA_seq = no_TGA_seq + "TGG"
        else:
            no_TGA_seq = no_TGA_seq + seq[(i*3):(i*3)+3]

    problem = DnaOptimizationProblem(
        sequence=str(no_TGA_seq),
        constraints=[AvoidPattern("BsaI_site"), AvoidPattern("NotI_site"), AvoidPattern("BbsI_site"), AvoidPattern("SapI_site"),
                     AvoidPattern(str(Seq.Seq(BsaI.site).reverse_complement())), AvoidPattern(str(Seq.Seq(NotI.site).reverse_complement())),
                     AvoidPattern(str(Seq.Seq(BbsI.site).reverse_complement())), AvoidPattern(str(Seq.Seq(SapI.site).reverse_complement())),
                     AvoidPattern(HomopolymerPattern("A", 7)),
                     AvoidPattern(HomopolymerPattern("T", 7)),
                     AvoidPattern(HomopolymerPattern("C", 7)),
                     AvoidPattern(HomopolymerPattern("G", 7)),
#                      AvoidPattern(RepeatedKmerPattern(2,5)),
                     EnforceTranslation()],
        objectives=[CodonOptimize(codon_usage_table=target_codon_usage, method='match_codon_usage',original_codon_usage_table=initial_codon_usage)])
    problem.resolve_constraints()
    problem.optimize()
    
    #To remove the sdterr (really annoying)
    sys.stderr = orig_stderr
    f.close()
    return str(problem.sequence)        
    
    
    
    
def create_codon_table(all_CDSs, number_table):
    
    """
    This function create a codon table.  
    
    ARGS:
    list of CDSs (list of strings)
    
    RETURN:
    A dictionnary in this format {'TTT': {'AA': 'F', 'AA_freq': 0.81, 'count': 10213}, 'TTC': {'AA': 'F', 'AA_freq': 0.18, 'count': 2375}...}  
    
    """
    
    import pandas as pd
    from Bio.Seq import Seq

    possible_nt = ['A','T','C','G']
    all_codons = [i+j+k for i in possible_nt for j in possible_nt for k in possible_nt]

    #Initialize the codon_table
    codon_table = {}
    for codon in all_codons:
        codon_table[str(codon)] = {
            'AA' : Seq(codon).translate(table=str(number_table))[0],
            'AA_freq' : None,
            'count' :0
        }

    #Construct the codon_table
    for cds in all_CDSs:
        for idx in range(int(len(cds)/3)):
            codon=cds[idx*3:idx*3+3]
            if codon in codon_table:
                codon_table[str(codon)]["count"] += 1                      
            else:
                break
                print('A big mistake is in this code. PLS verify the function create_codon_table in functions.py file')

    #make frequencies appear in dict:
    amino_acid_count = {}
    for codon in codon_table:
        if codon_table[codon]['AA'] not in amino_acid_count:
            amino_acid_count[codon_table[codon]['AA']] = codon_table[codon]['count']
        else:
            amino_acid_count[codon_table[codon]['AA']] += codon_table[codon]['count']


    #calculate frequencies:
    for codon in codon_table:
        codon_table[codon]['AA_freq'] = round(codon_table[codon]['count'] / amino_acid_count[codon_table[codon]['AA']],3)

    return codon_table  



def reverse_codon_table(codon_table):
    
    """
        This function reverse a codon table.
        
        ARGS: 
        codon table in this format: {'TTT': {'AA': 'F', 'AA_freq': 0.81, 'count': 10213}, 'TTC': {'AA': 'F', 'AA_freq': 0.18, 'count': 2375}...}
        
        RETURN:
        reverse codon table in this format: {'*': {'TAA': 0.81, 'TAG': 0.19}, 'F': {'TTT': 0.81, ...}...} 

    """

    rev_codon_usage = {codon_table[codon]['AA']:{} for codon in codon_table}

    for AA in rev_codon_usage:
        for codon in codon_table:
            if codon_table[codon]["AA"] == AA:
                rev_codon_usage[AA][codon] = codon_table[codon]["AA_freq"]
    return rev_codon_usage



#Other stuff

def count_codon(seq,codon):
    """
    Count the number of occurences of a specific codon
    
    Parameters:
        -seq (str) : sequence of DNA of length%3==0 (ex:'ATGCGGTAG')
        -codon (str): codon sequence in DNA (ex: 'ATG')
    
    Returns:
        - codon_count (int) : number of occurences of the codon in the sequence
    
    """
    
    assert len(seq)%3 == 0
    codon_count = 0
    for i in range(int(len(seq)/3)):
        if seq[(i*3):(i*3)+3] == codon:
            codon_count += 1
    return codon_count

def find_homopolymer(seq,length):
    """
    Find in a sequence contains an homopolymer of specified length
    
    Parameters:
        -seq (str) : sequence of DNA (ex:'ATGTTTTTTTTTTTTTTTTTAG')
        -length (int): homopolymer length to be searched (ex: 7)
    
    Returns:
        -True|False (Bool): return true if an homopolymer of specified length was found
    
    """
    for i in range(len(seq)-length+1):
        sub_seq = seq[i:i+length]
        sub_seq_uniq = set(sub_seq)
        if len(sub_seq_uniq) == 1:
            print("homopolymer found :",sub_seq,"position :",i)
            return True
    return False

def count_homopolymer(seq,length):
    """
    Count number of occurences of homopolymer of specified length in a given sequence
    
    Parameters:
        -seq (str) : sequence of DNA (ex:'ATGTTTTTTTTTTTTTTTTTAG')
        -length (int): homopolymer length to be searched (ex: 7)
    
    Returns:
        -count (int): number of occurences of homopolymer of specified length in a given sequence
    
    """
    count = 0
    for i in range(len(seq)-length+1):
        sub_seq = seq[i:i+length]
        sub_seq_uniq = set(sub_seq)
        if len(sub_seq_uniq) == 1:
            count += 1
    return count

def codon_pc_dif(cds1, cds2):
    """
    Compute the percentage of codons different between 2 sequences
    
    Parameters:
        -cds1 (str) : sequence of DNA of length%3==0 (ex:'ATGCGGTAG')
        -cds2 (str) : sequence of DNA of length%3==0 (ex:'ATGCGGTAG')
    
    Returns:
        -pc_dif_score (float): percentage of codons different between 2 sequences
    
    """
    assert len(cds1) == len(cds2)
    assert len(cds1)%3 == 0
    
    cds1 = str(cds1)
    cds2 = str(cds2)
    
    dif_score = 0
    
    for i in range(int(len(cds1)/3)):
        if cds1[(i*3):(i*3)+3] != cds2[(i*3):(i*3)+3]:
            dif_score += 1
    
    return (dif_score/(len(cds1)/3))*100
                
def find_right_coordinates(index,restricted_region,max_overlap_len,window_len,prev,current_left_border,verbose = False):
    """
    Calculates the position where a primer can be placed and pick a specific position. 
    
    Parameters:
        index: The higher the index is, smaller the length of the g_block is . If a primer cannot be found in a region, shift to another with index.. (int)
        restricted_region =  Contains all the region where a primer cannot be placed. (list)
        max_overlap_len =  The maximum overlap between two g-block for homologous recombinaison in yeast. (int > 0)
        window_len = The core length of the gblock. The g-block length - homology. (int)
        prev =  the position of the previous right primer. (int)
        current_left_border = See the jupyter notebook for more info. Generally it is left_coordinates[0]. (int)
        
    
    Returns:
        A list of coordinates ([start, end]). It is the coordinates for the right end of the g_block. 
    
    """
    if prev == 0:
        if verbose:print (str(index)+" fois!")

        #trouver la region qui pourrait etre celle d'interet: ajouter une region interdite si prev trouve apres 
        #Notes: il ne devrait pas y avoir de regions interdites de plus de 1,8kb: un message d'erreur par possible_job est envoye a ce moment-la. 
        #(POSSIBLE_JOB N'EST PAS ACHEVE)
        little_restriction = [element for element in restricted_region if int(element.split(':')[1]) > 0 and 
                              int(element.split(':')[0]) < 0+window_len+max_overlap_len]

        #lister toutes les places qui peuvent servir a primer3 pour trouver les primers possibles de droite: 
        liste_coordinates_sequence_template = []
        
        if len(little_restriction) == 0:
            liste_coordinates_sequence_template.append([window_len,max_overlap_len+window_len])
        
        else:
            for idx, element in enumerate(little_restriction):
                right = int(element.split(':')[1])
                left = int(element.split(':')[0])

                #if the first restricted region starts after max_overlap_len, [0:left_restricted_region_border] can be a gblock.
                if idx == 0 and left >= max_overlap_len:
                    liste_coordinates_sequence_template.append([(prev-current_left_border)+current_left_border+max_overlap_len,left])

                #for all the other restricted region between 0 to 1800:
                if idx < len(little_restriction)-1:
                    next_left = int(little_restriction[idx+1].split(':')[0])
                    if next_left - right >= max_overlap_len:
                        add = [right,next_left]
                        liste_coordinates_sequence_template.append(add)

                #if the right border of the last restricted region is >200 bp long, [right_border, max_overlap_len+window_len+0] can be a gblock.
                if idx == len(little_restriction)-1:
                    if (max_overlap_len + window_len ) - right >= max_overlap_len:
                        liste_coordinates_sequence_template.append([right,max_overlap_len+window_len])

            for coordinates in liste_coordinates_sequence_template:
                if coordinates[1] - coordinates[0] < max_overlap_len :
                    liste_coordinates_sequence_template.remove(coordinates)
                    
        if verbose:print('je suis toutes les zones pouvant accueillir des primers: '+ str(liste_coordinates_sequence_template))
        all_possible_places_for_primer = generate_all_primer_places(liste_coordinates_sequence_template, max_overlap_len)
        if verbose:print('je suis tous les places possibles pour accueillir des primers: '+ str(all_possible_places_for_primer))
        coordinates_sequence_template = all_possible_places_for_primer[-(index)]
            
    #In the general case:
    if prev != 0:
        if verbose:print (str(index)+" fois!")
        
        #lister toutes les places qui peuvent servir a primer3 pour trouver les primers possibles de droite:
        little_restriction = [element for element in restricted_region if int(element.split(':')[1]) > current_left_border and 
                              int(element.split(':')[0]) < current_left_border+window_len+max_overlap_len]
        
        if verbose:print('je suis little_restriction: '+str(little_restriction))
        
        right_coordinates_sequence_template = []
        
        if len(little_restriction) == 0:
            right_coordinates_sequence_template.append([current_left_border+window_len,max_overlap_len+current_left_border+window_len])
        
        else:
            for idx, element in enumerate(little_restriction):

                right = int(element.split(':')[1])
                left = int(element.split(':')[0])

                if idx == 0 and left >= max_overlap_len + current_left_border:
                    #right_coordinates_sequence_template.append([current_left_border+max_overlap_len,left])
                    right_coordinates_sequence_template.append([(prev-current_left_border)+current_left_border+max_overlap_len,left])
                if idx < len(little_restriction)-1:
                    next_left = int(little_restriction[idx+1].split(':')[0])
                    if next_left - right >= max_overlap_len:
                        add = [right,next_left]
                        right_coordinates_sequence_template.append(add)
                if idx == len(little_restriction)-1:
                    #if the right border of the last restricted region + 200 is in the limit of the gblock:
                    if (max_overlap_len + window_len + current_left_border) - right >= max_overlap_len:
                        #right_coordinates_sequence_template.append([window_len+current_left_border,window_len+max_overlap_len+current_left_border])
                        right_coordinates_sequence_template.append([right,window_len+max_overlap_len+current_left_border])
                        
            for coordinates in right_coordinates_sequence_template:
                if coordinates[1] - coordinates[0] < max_overlap_len :
                    right_coordinates_sequence_template.remove(coordinates)
                    
        if verbose:print("je suis right_coordinates_sequence_template : "+str(right_coordinates_sequence_template))
        all_possible_places_for_primer = generate_all_primer_places(right_coordinates_sequence_template, max_overlap_len)
        if verbose:print('je suis tous les places possibles pour accueillir des primers: '+ str(all_possible_places_for_primer))
        if verbose:print('je suis index:',index)    
        coordinates_sequence_template = all_possible_places_for_primer[-(index)]
            
    return coordinates_sequence_template


def generate_all_primer_places(right_coordinates_sequence_template, max_overlap_len):
    """
    Creates all possible 200bp that can be used as right sequence template. 
    
    Parameter:
        right_coordinates_sequence_template: a list of coordinates that defines the limits of where a primer can bind (based on restricted regions). 
    
    Returns:
        a list of all the 200bp (max overlap) possible based on the interval provide in input. 
    
    """
    ordered_list = []
    for element in right_coordinates_sequence_template:
        new_list = []
        decrease_number = element[1]
        while decrease_number - element[0] >= 2*max_overlap_len:
            new_list.insert(0,[int(decrease_number - max_overlap_len),int(decrease_number)])
            decrease_number = decrease_number - max_overlap_len
        new_list.insert(0,[int(element[0]),int(element[1])])
        ordered_list.insert(1,new_list)

    output_list = sorted([sub_element for element in ordered_list for sub_element in element], key=lambda x: x[0]) 

    return output_list