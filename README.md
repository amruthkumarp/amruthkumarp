makeblastdb -in HSA_Proteins.fasta -dbtype prot -out HSA_Proteins_DB
blastp -query ThiQ.fasta -db HSA_Proteins_DB -out ThiQ_vs_HSA.blast -evalue 1E-04 -outfmt 6
Query ID
Subject ID
% Identity
Alignment Length
Mismatches
Gap Opens
Query Start
Query End
Subject Start
Subject End
E-value
Bit Score
cut -f2 ThiQ_vs_HSA.blast | sort | uniq -c | awk '$1 > 1'
makeblastdb -in HSA_Proteins.fasta -dbtype prot -out HSA_Proteins_DB
psiblast -query SRY.fasta -db HSA_Proteins_DB -num_iterations 3 -evalue 1E-04 -out SRY_PSIBLAST_iter3.txt -outfmt 6 -num_threads 4 -save_pssm_after_last_round
grep -A 1000 "Results from round 3" SRY_PSIBLAST_iter3.txt | tail -n +2 > SRY_paralogs_iter3.txt
cat SRY_paralogs_iter3.txt
blastp -query ThiQ.fasta -db HSA_Proteins_DB -out ThiQ_vs_HSA.blast -evalue 1E-04 -outfmt 6 -max_target_seqs 1 -num_threads 4
needle -asequence ThiQ.fasta -bsequence Human_ThiQ.fasta -gapopen 10 -gapextend 0.5 -outfile Global_Alignment.txt
water -asequence ThiQ.fasta -bsequence Human_ThiQ.fasta -gapopen 10 -gapextend 0.5 -outfile Local_Alignment.txt
cat Global_Alignment.txt
cat Local_Alignment.txt
blastp -query SRY.fasta -db nr -out SRY_BLAST_results.txt -evalue 1E-05 -outfmt 6 -max_target_seqs 1000 -num_threads 4
awk '$3 >= 60 && $3 < 70' SRY_BLAST_results.txt | head -n 5 > SRY_60_69.txt
awk '$3 >= 70 && $3 < 80' SRY_BLAST_results.txt | head -n 5 > SRY_70_79.txt
awk '$3 >= 80 && $3 < 90' SRY_BLAST_results.txt | head -n 5 > SRY_80_89.txt
awk '$3 >= 90 && $3 <= 100' SRY_BLAST_results.txt | head -n 5 > SRY_90_100.txt
cut -f2 SRY_60_69.txt > IDs_60_69.txt
cut -f2 SRY_70_79.txt > IDs_70_79.txt
cut -f2 SRY_80_89.txt > IDs_80_89.txt
cut -f2 SRY_90_100.txt > IDs_90_100.txt

clustalo -i SRY_All_Homologues.fasta -o SRY_MSA.aln --outfmt=clu --threads=4
cat SRY_MSA.aln
cat SRY_60_69.fasta SRY_70_79.fasta SRY_80_89.fasta SRY_90_100.fasta > SRY_All_Homologues.fasta
tblastn -query SRY_All_Homologues.fasta -db nt -out SRY_DNA_Hits.txt -evalue 1E-05 -outfmt 6 -max_target_seqs 1 -num_threads 4
cut -f2 SRY_DNA_Hits.txt > DNA_Accession_IDs.txt
blastdbcmd -db nt -entry_batch DNA_Accession_IDs.txt -out SRY_Homologues_DNA.fasta
cat SRY_Homologues_DNA.fasta
grep -c ">" SRY_Homologues_DNA.fasta
cat SRY_Homologues_DNA.fasta > SRY_Motif_Input.fasta
clustalo -i SRY_All_Homologues.fasta -o SRY_MSA.aln --outfmt=clu --threads=4
awk '($1~"^ATOM$" && $5~"^B$")' 1w9p.pdb  > 1w9p_b.pdb
awk '($1~"^ATOM$" && $5~"^B$")' 3aqu.pdb  > 3aqu_b.pdb
awk '($1~"^ATOM$" && $5~"^B$")' 1i54.pdb  > 1i54_b.pdb
mustang -i 1w9p_b.pdb 1i54_b.pdb -o newres








Hmm
ifconfig | grep 192 > ip.txt
date > date.txt
Copy the command below and press ENTER
export HISTTIMEFORMAT="%F %T %z "
history 50 > history.txt
last > last.txt

before doing this do multiple sequence alignment
clustalo 
clustalo -i msa_input.fasta -o example_alignment.sto --outfmt=st

hmmbuild	Builds a profile HMM from a multiple sequence alignment
hmmsearch	Searches a profile HMM against a sequence database
hmmscan	Searches sequences against an HMM database
hmmpress	Prepares an HMM database for scanning
hmmalign	Aligns sequences to a profile HMM

hmmbuild
hmmbuild my_model.hmm aligned_sequences.sto

hmmsearch
hmmsearch my_model.hmm proteins.fasta > results.txt
hmmsearch --tblout search_table.txt example.hmm proteins.fasta


hmmscan
hmmscan Pfam-A.hmm sequences.fasta > scan_results.txt

hmmpress
hmmpress Pfam-A.hmm(before feeding it to hmmscan)

hmmalign
hmmalign my_model.hmm unaligned.fasta > aligned.sto
hmmalign -o output.sto profile.hmm sequences.fasta


hmmstat Pfam-A.hmm > Pfam_profile.txt
hmmfetch Pfam-A.hmm PF03790.16 > retrieved_profile.hmm
hmmpress Pfam-A.hmm
hmmscan Pfam-A.hmm P24740.faa > P24740_in_Pfam

phmmer
phmmer query.fasta target_proteins.fasta > output.txt

jackhmmer

jackhmmer <filname> <database> > <outputfile> 
