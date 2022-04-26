#!/usr/bin/env perl

use strict;
use warnings;
#use diagnostics;
use LWP::UserAgent;
use Scalar::Util qw(openhandle);
use Bio::SeqIO;
use Bio::Seq;
use Data::Dumper;
use IO::File;
use SWISS::Entry;
use File::Basename;
use List::Util qw(min max sum sum0);
use List::MoreUtils qw(uniq);
use Scalar::Util qw(openhandle);

# !!! THIS SCRIPT HAS NOT BEEN RUN"AS SUCH". IT IS A COLLECTION OF ALL BITS AND PIECES USED TO CREATE THE FINAL DATABASE




# https://www.biostars.org/p/238376/

#get ko00000.keg file for all K numbers Dowload text ko00000.keg file from https://www.genome.jp/kegg-bin/get_htext?ko00000.keg
#https://www.genome.jp/kegg-bin/download_htext?htext=ko00000.keg&format=htext&filedir=

#get orthology_uniprot.list from linkDB kegg

### up:O50214 Styrene monooxygenase StyA (EC:1.14.14.11) ko:14481 was not present in "orthology_uniprot.list" 
#  for some reason, causing the script to stop. So I added it manualy to the file  with [echo 'ko:K14481     up:O50214   reverse' >> orthology_uniprot.list]
###  [ko:K18541    up:Q84DC4       reverse]     mandelamide amidase [EC:3.5.1.86]
###  [ko:K01628    up:A0A1H6VLF3   reverse]	L-fuculose 1-phosphate aldolase
###  [ko:K00058    up:A0A1H6U3H5   reverse]   	D-3-phosphoglycerate dehydrogenase / 2-oxoglutarate reductase [EC:1.1.1.95 1.1.1.399]
###  [ko:K00058    up:A0A1H6VQX1   reverse]  	D-3-phosphoglycerate dehydrogenase / 2-oxoglutarate reductase [EC:1.1.1.95 1.1.1.399]
###  [ko:K00757    up:A0A1H6VWW5   reverse]  	uridine phosphorylase [EC:2.4.2.3]
###  [ko:K00768    up:A0A1H6RC40   reverse]  	nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase [EC:2.4.2.21]
###  [ko:K00945    up:A0A1H6W7I6   reverse]  	CMP/dCMP kinase [EC:2.7.4.25]
###  [ko:K00973    up:A0A1H6XRP5   reverse]  	glucose-1-phosphate thymidylyltransferase [EC:2.7.7.24]
###  [ko:K01101    up:A0A1H6RGM0   reverse]  	4-nitrophenyl phosphatase [EC:3.1.3.41]
###  [ko:K01486    up:A0A1H6RWE5   reverse]  	adenine deaminase [EC:3.5.4.2]
###  [ko:K01628    up:A0A1H6WAW7   reverse] 	L-fuculose-phosphate aldolase [EC:4.1.2.17]
###  [ko:K02304    up:A0A1H6XB12   reverse]  	precorrin-2 dehydrogenase / sirohydrochlorin ferrochelatase [EC:1.3.1.76 4.99.1.4]
###  [ko:K03752    up:A0A1H6RNW2   reverse]  	molybdenum cofactor guanylyltransferase [EC:2.7.7.77]
###  [ko:K03799    up:A0A1H6SPC0   reverse]  	heat shock protein HtpX [EC:3.4.24.-]
###  [ko:K03841    up:A0A1H6V201   reverse]  	fructose-1,6-bisphosphatase I [EC:3.1.3.11]
###  [ko:K05520    up:A0A1H6W4C7   reverse]  	protease I [EC:3.5.1.124]
###  [ko:K05577    up:P0DPI1   reverse] 	NAD(P)H-quinone oxidoreductase subunit 5 [EC:7.1.1.2]
###  [ko:K06162    up:A0A1H6XE13   reverse]  	alpha-D-ribose 1-methylphosphonate 5-triphosphate diphosphatase [EC:3.6.1.63]
###  [ko:K06163    up:A0A1H6X1Z1   reverse] 	alpha-D-ribose 1-methylphosphonate 5-phosphate C-P lyase [EC:4.7.1.1]
###  [ko:K06164    up:A0A1H6X213   reverse]  	alpha-D-ribose 1-methylphosphonate 5-triphosphate synthase subunit PhnI [EC:2.7.8.37]
###  [ko:K11780    up:A0A1H6RFK5   reverse]  	7,8-didemethyl-8-hydroxy-5-deazariboflavin synthase [EC:4.3.1.32]
###  [ko:K13280    up:A0A1H6QTY3   reverse]  	signal peptidase I [EC:3.4.21.89]

### K15391  CPAA; cyclopiazonic acid synthetase, hybrid polyketide synthase / nonribosomal peptide synthetase 
#  did not have an entry in Uniprot. removed this manually from ko00000.keg

### K20926  cap3B cellobiuronic acid synthase [EC:2.4.1.-] removed this manually from ko00000.keg
### K22634  salviol synthase [EC:1.14.14.62] CYP71BE52 removed this manually from ko00000.keg
### K22635  CYP76AH3; ferruginol monooxygenase / sugiol synthase removed this manually from ko00000.keg
### K22636  CYP76AH24; ferruginol synthase / ferruginol monooxygenase removed this manually from ko00000.keg
### K22637  CYP76AK1; 11-hydroxysugiol 20-monooxygenase removed this manually from ko00000.keg
### K22638  CYP76AK6; carnosic acid synthase removed this manually from ko00000.keg
### K22639  CYP92C5; trimethyltridecatetraene/dimethylnonatriene synthase removed this manually from ko00000.keg
### K22640  CYP92C6; trimethyltridecatetraene synthase removed this manually from ko00000.keg

### The following samples were not part of the metabolic section of ko00000.keg, but did appear in the ortho_uniprot list.
# Added those to META.KO.txt [echo 'D      K00344  qor, CRYZ; NADPH2:quinone reductase' >> META.KO.txt]
### K00344 D      K00344  qor, CRYZ; NADPH2:quinone reductase
### K03671 D      K03671  trxA; thioredoxin 1
### K03676 D      K03676  grxC, GLRX, GLRX2; glutaredoxin 3
### K11991 D      K11991  tadA; tRNA(adenine34) deaminase
### K03439 D      K03439  trmB, METTL1; tRNA (guanine-N7-)-methyltransferase
### K08483 D      K08483  PTS-EI.PTSI, ptsI; phosphotransferase system, enzyme I, PtsI
### K11183 D      K11183  PTS-HPR.FRUB, fruB, fpr; phosphocarrier protein FPr
### K07219 D      K07219  K07219; putative molybdopterin biosynthesis protein
### K03741 D      K03741  ARSC2, arsC; arsenate reductase
### K04015 D      K04015  NrfD; protein NrfD
### K00153 D      K00153  E1.1.1.306; S-(hydroxymethyl)mycothiol dehydrogenase [EC:1.1.1.306]
### K03518 D      K03518  coxS; aerobic carbon-monoxide dehydrogenase small subunit
### K03519 D      K03519  coxM, cutM; aerobic carbon-monoxide dehydrogenase medium subunit
### K03520 D      K03520  coxL, cutL; aerobic carbon-monoxide dehydrogenase large subunit
### K08356 D      K08356  aoxB; arsenite oxidase large subunit
### K08355 D      K08355  aoxA; arsenite oxidase small subunit
### K17050 D      K17050  clrA, serA; complex iron-sulfur molybdoenzyme family reductase subunit alpha
### K17051 D      K17051  clrb, serB; complex iron-sulfur molybdoenzyme family reductase subunit beta
### K17052 D      K17052  clrC, serC; complex iron-sulfur molybdoenzyme family reductase subunit gamma
### K09162 D      K09162  cld; chlorite dismutase
### K01539 D      K01539  ATP1A; sodium/potassium-transporting ATPase subunit alpha
### K01539 D      K02945  RP-S1, rpsA; small subunit ribosomal protein S1
### K01539 D      K02487  pilL; type IV pili sensor histidine kinase and response regulator
### K01539 D      K11005  hlyA; hemolysin A
### K01539 D      K01090  E3.1.3.16; protein phosphatase




### The following samples were added to META.KO.txt because they form a cluster that fits metabolically into the represented K number
### but they do not have a K number yet or are not recognised by Kegg as such. They are, however, different enough from the canonical homolog
### to not be included as the canonical homolog.
 
# The metadata in the final HMM was changed to-> K03385  otr; octaheme c-type cytochrome, tetrathionate reductase family
### K03385  otr; octaheme c-type cytochrome, tetrathionate reductase family => nrfA like.



#my $KO_file="/scratch2/databases/metascan/ko00000.keg"; NOT CLEANED YET^^^^^^^^^^^^^^^^^^^^^^
my $KO_file="ko00000.keg";
my $ortho_uniprot="orthology_uniprot.list";

# select metabolic genes from KO

print "Selecting metabolic genes from KO\n";
system("sed -n '/A09100/,/A09120/p' $KO_file > META.KO.txt");

#grep  KO numbers and create a hash with Knumbers and gene names

print "Counting KO with gene name from metabolic section of Kegg list\n";
open my $stats, '>', "stats.txt" or die $!;
my %keg;
my $count=0;
open (KO, "<META.KO.txt") or die "Can't open META.KO.txt";
while (my $line = <KO>){
   chomp $line;
   if ($line =~ /^D/) {
      my ($D, $ko, $desc) = split(/\s+/, $line, 3);
      $keg{$ko}=$desc;
      $count++;
   }
}
print "->$count\n";
print {$stats} "Counting KO with gene name from metabolic section of Kegg list\n->$count\n";

close (KO);

#create a hash with Knumbers and uniprot numbers

print "Creating a hash with KO and Uniprot numbers\n";


my %ortho_rev;
open (OU, "<$ortho_uniprot") or die "Can't open $ortho_uniprot";
while (my $line = <OU>){
   chomp $line;
   my ($D, $ko,$a, $uniprot) = split(/[:\s]+/, $line, 5);
   push (@{$ortho_rev{$uniprot}}, $ko);
}
close (OU);

my %ortho;
open (OU, "<$ortho_uniprot") or die "Can't open $ortho_uniprot";
while (my $line = <OU>){
   chomp $line;
   my ($D, $ko,$a, $uniprot) = split(/[:\s]+/, $line, 5);
   push (@{$ortho{$ko}}, $uniprot);
}
close (OU);


#create file with KO numbers and gene names

print "Creating a file with KO and gene names:Kegg_name.txt\n";

open my $Kn , '>', "Kegg_name.txt" or die $!;
my $count2=0;
for my $key (sort keys %keg) {
   print {$Kn} "$key\t$keg{$key}\n"; $count2++;
}

print "->$count2\n"; 
print {$stats} "Creating a file with KO and gene names:Kegg_name.txt\n->$count2\n";
#create files with uniprot numbers for each metabolic K numbers

print "Creating files with uniprot numbers for each metabolic K numbers\n";

mkdir "up_3";
mkdir "up_3_embl";
mkdir "up_3_fasta";
mkdir "up_012";
mkdir "up_3_trim";

open my $under_3, '>', "under3.txt" or die $!;

foreach my $K (sort keys %keg){
   if (scalar @{$ortho{$K}} >= 3){
      open my $KOU, '>', "up_3/$K" or die $!;
      print {$KOU} "@{$ortho{$K}}\n";
   }
   if (scalar (@{$ortho{$K}}) < 3){
      open my $KOU, '>', "up_012/$K" or die $!;
      print {$KOU} "@{$ortho{$K}}\n";
   }
   if (scalar (@{$ortho{$K}}) < 3)  {
      print {$under_3} "$K\t$keg{$K}; Not enough entries for an HMM profile\n"
   }
}

#download embl files from uniprot

print "Downloading TrEMBL-prot files from EMBL\n";

opendir (DH, "up_3/") or die "Can't opendir up_3/ for reading '$!'";
my @files =  readdir(DH);
my @sorted_files = sort { $a cmp $b } @files;
closedir (DH);

my $count3=0;
foreach my $file (@sorted_files){
   if($file eq "." || $file eq ".." || -e "up_3_embl/$file.embl"){
      next;
   }
   $count3++;
   print "$file\t$count3\n";

   open my $embl_fh , '>', "up_3_embl/$file.embl"; 
   my $list = "up_3/$file"; # File containing list of UniProt identifiers.

   my $base = 'http://www.uniprot.org';
   my $tool = 'uploadlists';

   my $contact = 'gcremers@science.ru.nl'; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
   my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact",
                                   keep_alive => 1);

   push @{$agent->requests_redirectable}, 'POST';
   my $response = $agent->post("$base/$tool/",
                            [ 'file' => [$list],
                              'format' => 'txt',
                              'from' => 'ACC+ID',
                              'to' => 'ACC',
                            ],
                              'Content_Type' => 'form-data');
   while (my $wait = $response->header('Retry-After')) {
      print STDERR "Waiting ($wait)...\n";
      sleep $wait;
      $response = $agent->get($response->base);
   }

   $response->is_success ?
   print {$embl_fh} $response->content :
   die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
}

# convert files to fasta
print "Converting to fasta\n";

opendir (EMBL, "up_3_embl/") or die "Could not open up_3_embl/ for reading '$!'\n";
my @embls = grep {$_ =~ /\.embl$/} readdir (EMBL);
my @sorted_embls = sort { $a cmp $b } @embls;

foreach my $embl (@sorted_embls){
   if(-e "up_3_fasta/$embl.fasta"){
      next;
   }

   my $fh = new IO::File "up_3_embl/$embl" or 
    die "Cannot open input file up_3_embl/$embl: $!";
   open my $fsa_fh, '>', "up_3_fasta/$embl.fasta" or die $!;


   $/ = "\n\/\/";
   while(<$fh>) {
      s/\r//g;
      (my $entry_txt = $_) =~ s/^\s+//;
      next unless $entry_txt;
      $entry_txt .= "\n";
      my $entry = SWISS::Entry->fromText(e $entry_txt );
      print {$fsa_fh} $entry->toFasta();
   }
}
closedir (EMBL);

# Because Bioperl does not retrieve the ID from these samples/files and SWISS::KNIFE not the length,
# trimming the fastas on length has to be done afterwards :/ Besides, there is something in 
# some embl files that causes the code to crash, that SWISS::KNIFE does not seem to matter.   


print "Trimming fasta\n";


opendir (TRIM, "up_3_fasta/") or die "Could not open up_3_fasta/ for reading '$!'\n";
my @trims = grep {$_ =~ /\.fasta$/} readdir (TRIM);
my @sorted_trims = sort { $a cmp $b } @trims;

foreach my $trim (@sorted_trims){
   my $base = basename("$trim",  ".embl.fasta");
   if(-e "up_3_trim/$base.trim.fasta"){
      next;
   }
   my $seq_in = Bio::SeqIO->new( -format => 'fasta',
                                 -file   => "up_3_fasta/$base.embl.fasta",
                            );
   my $seq_out = Bio::SeqIO->new( -file   => ">up_3_trim/$base.trim.fasta",
                                  -format => 'fasta',
                             );
   open my $trimstat_fh, '>>', "trimstatistics.txt" or die $!;

   # loads the whole file into memory - be careful
   # if this is a big file, then this script will
   # use a lot of memory
   my $startcount=0;
   my @seq_array;
   my @check_array;
   while ( my $seq = $seq_in->next_seq() ) {
      $startcount++;
      my $doublecheck=$seq->{'primary_seq'}->{'seq'};
      if ( grep (/^$doublecheck$/, @check_array) ) { #dereplication
         next; 
      }
      else {
         push (@check_array, $seq->{'primary_seq'}->{'seq'}); 
         push (@seq_array,$seq);
      }   
   }
   # now do something with these. First sort by length,
   # find the average and median lengths and print them out

   @seq_array = sort { $a->length <=> $b->length } @seq_array;

   my $total = 0; 
   my $count = 0; 
   for my $seq ( @seq_array ) { 
      $total += $seq->length;
      $count++;
   }
 
   my $upper_limit= $total/$count*1.5; #set upper size limit to 150% of mean
   my $lower_limit= $total/$count*0.6; #set lower limit to 60% length of mean
   
   print {$trimstat_fh} "$base\n"; 
   print {$trimstat_fh} "Total protein: ", $startcount, "\t", "total AA: $total\n";
   print {$trimstat_fh} "Unique sequences:  ", scalar @seq_array, "\n";
   print {$trimstat_fh} "Mean length: ", $total/$count, " Median: ", 
         $seq_array[$count/2]->length, "\n";
   print  {$trimstat_fh} "Upper limit: $upper_limit\t", "Lower limit: $lower_limit\n";
 

   my $counttrim=0; 
   for my $seq ( @seq_array ) {
      if ($seq->length > $lower_limit and $seq->length < $upper_limit){
          $counttrim ++;
      }
   }
   if ($counttrim < 3)  {
      print {$under_3} "$base\t$keg{$base}: after trimming and dereplication. See trimstatistics.txt for more info.\n";
   }   

   for my $seq (@seq_array) {
      if ($counttrim >= 3) { #if at least 3 sequence are left after trimming
         if ($seq->length > $lower_limit and $seq->length < $upper_limit){
             $seq_out->write_seq($seq);
         }
      }
      if ($counttrim < 3){ #otherwise use untrimmed set
         $seq_out->write_seq($seq);
      }
   }
   if ($counttrim < 3 and scalar @seq_array >= 3){ 
         print {$trimstat_fh} "Protein left after trimming: $counttrim;-> lower than 3: use untrimmed set\n";
   }
   if ($counttrim >= 3){ 
         print {$trimstat_fh} "Protein left after trimming: $counttrim\n";
   }
   if ($count < 3){
         print {$trimstat_fh} "Unique Proteins left: $count; Not enough unique sequences for an Hmm profile:-> Discard set\n";
   }
   print {$trimstat_fh} "<----->\n";
} 
closedir (TRIM);

print "Clustering\n";
mkdir "up_3_clusters";

# Some Knumbers (like K07540 and K05301 contain also subunits) These should cluster 
# together during clustering resulting in different HMM-profiles. But some (like K05301 
# are all too different to cluster and will result in no clusters and thus no HMM-profile.

# K10966 has three identical sequences. Which apparently causes Linclust to crash. dereplicate in trimming

open my $no_hhm_fh, '>>', "HMM_empty.txt" or die $!;

opendir (CLUST, "up_3_trim/") or die "Could not open up_3_trim/ for reading '$!'\n";
my @clusts = grep {$_ =~ /\.fasta$/} readdir (CLUST);
my @sorted_clust = sort { $a cmp $b } @clusts;

foreach my $clust (@sorted_clust){
   my $base = basename("$clust",  ".trim.fasta");
   my $KMER="160";
   if( -d -e "up_3_clusters/dir.$base"){
      next;
   }

   my $ID="0.5";
 
   mkdir "up_3_clusters/dir.$base";
   mkdir "up_3_clusters/dir.$base/clusters";
   mkdir "up_3_clusters/dir.$base/clusters/align";
   mkdir "up_3_clusters/dir.$base/clusters/hmm";

   system("/scratch/MMseqs2/build/bin/mmseqs createdb -v 0 up_3_trim/$base.trim.fasta $base.DB");
   system("/scratch/MMseqs2/build/bin/mmseqs linclust -v 0 --kmer-per-seq $KMER --min-seq-id $ID --similarity-type 1 --sub-mat /scratch/MMseqs2/data/blosum80.out --cluster-mode 2 --cov-mode 0 -c 0.7 $base.DB clu tmp.$base.$KMER ");
   system("/scratch/MMseqs2/build/bin/mmseqs createseqfiledb -v 0 $base.DB clu clu_seq");
   system("/scratch/MMseqs2/build/bin/mmseqs result2flat -v 0 $base.DB $base.DB clu_seq clu_seq.fasta");
   system("/scratch/MMseqs2/build/bin/mmseqs createtsv -v 0 $base.DB $base.DB clu clu.tsv");
   #system("/scratch/MMseqs2/build/bin/mmseqs result2repseq $base.DB clu DB_clu_rep");
   #/scratch/MMseqs2/build/bin/mmseqs result2flat DB DB DB_clu_rep DB_clu_rep.fasta  --use-fasta-header && \
   system("mv c* up_3_clusters/dir.$base && mv $base.DB* up_3_clusters/dir.$base"); # mv D* up_3_clusters/dir.$base.$KMER &&
   system("rm -rf tmp.$base.$KMER");

   my $tsv="up_3_clusters/dir.$base/clu.tsv";
   my %tsv;
   open (TSV, "<$tsv") or die "Can't open $tsv";
      while (my $line = <TSV>){
         chomp $line;
         my ($bin, $up) = split(/\t/, $line, 2);
         push(@{$tsv{$bin}}, $up);
      }
   close (TSV);

   my @miscs;
   foreach my $bin (keys %tsv) {
      if (scalar @{$tsv{$bin}} < 3){ 
         push(@miscs, $bin);
      }
   }
   
   my @sorted_miscs = sort { $a cmp $b } @miscs;
   my $misc= $sorted_miscs[0];
   
   my $count=0;
   foreach my $bin (keys %tsv) {
      if (scalar @{$tsv{$bin}} >= 3){
         $count++;
         my $seqin = Bio::SeqIO->new(-file=>"<up_3_clusters/dir.$base/clu_seq.fasta", -format=>'fasta');
         my $seqout = Bio::SeqIO->new(-file=>">up_3_clusters/dir.$base/clusters/$bin.$count.faa", -format=>'fasta');
         while ( my $seq = $seqin->next_seq() ) {
            foreach my $id (@{$tsv{$bin}}) {
               if ($seq->display_id =~ /\|$id/i) {
                  $seqout->write_seq($seq);
               }
            }
         }
      }
      if (scalar @{$tsv{$bin}} < 3){ #gobble up all unused sequences into 1 big file to keep the script going and not lose to much data
         my $seqin = Bio::SeqIO->new(-file=>"<up_3_clusters/dir.$base/clu_seq.fasta", -format=>'fasta');
         my $seqout = Bio::SeqIO->new(-file=>">>up_3_clusters/dir.$base/clusters/$misc.misc.faa", -format=>'fasta'); #one binname only.
         while ( my $seq = $seqin->next_seq() ) {
            foreach my $id (@{$tsv{$bin}}) { 
               if ($seq->display_id =~ /\|$id/i) {
                  $seqout->write_seq($seq);
               }
            }
         }
      }
   }
 
   # alignment of fasta files
   print "Alignment of fasta files\n";

   opendir (FSA, "up_3_clusters/dir.$base/clusters/") or die "Could not open up_3_clusters/dir.$base/clusters/for reading '$!'\n";
   my @fastas = grep {$_ =~ /\.faa$/} readdir (FSA);
   my @sorted_fastas = sort { $a cmp $b } @fastas;
  
   foreach my $fasta (@sorted_fastas){
      if(-e "up_3_clusters/dir.$base/clusters/align/$fasta.msa"){
         next;
      }
      my $count=0;
      open FASTA, "up_3_clusters/dir.$base/clusters/$fasta" or die "can not open fasta";
      while (my $line = <FASTA>){
         chomp $line;       
         if ($line =~ /^>/) {
            $count++;
        }
      }
      close FASTA;
      if ($count >= 3){ #HERE IS THE ACTUAL TRIMMING ON >3 unique sequences
         system("mafft --quiet --anysymbol up_3_clusters/dir.$base/clusters/$fasta > up_3_clusters/dir.$base/clusters/align/$fasta.msa");
      }   
   }
   closedir (FSA);

   print "Creating HMM profiles\n";
   my $hmmcount=0;
   opendir (HMM, "up_3_clusters/dir.$base/clusters/align") or die "Could not open up_3_clusters/dir.$base/clusters/align for reading '$!'\n";
   my @alis = grep {$_ =~ /\.msa$/} readdir (HMM);
   my @sorted_alis = sort { $a cmp $b } @alis;

   foreach my $ali (@sorted_alis){
      if(-e "up_3_clusters/dir.$base/clusters/hmm/$ali.hmm"){
         next;
      }
      if(-e "up_3_clusters/dir.$base/clusters/align/$ali"){
         system("hmmbuild up_3_clusters/dir.$base/clusters/hmm/$ali.hmm up_3_clusters/dir.$base/clusters/align/$ali");
         $hmmcount++;
      }
   }

   if ($hmmcount == 0 ){
      print {$no_hhm_fh} "$clust\n";
   }
   if ($hmmcount > 0 ){
      system("cat up_3_clusters/dir.$base/clusters/hmm/*.hmm > up_3_clusters/dir.$base/clusters/hmm/$base.hmm ");
      system("cat up_3_clusters/dir.$base/clusters/hmm/$base.hmm >> meta.hmm");
   }
   closedir (HMM);
}
closedir (CLUST);

my $meta="meta.hmm";

open my $in,  '<',  $meta      or die "Can't read old file: $!";
open my $out, '>', "$meta.new" or die "Can't write new file: $!";
#print {$out} "#CYCLE METAB\n\n";
foreach  (<$in>)  {
   if ($_ =~ /^NAME/) {
      chomp $_;
      my ($name, $fname) = split(/\s+/, $_, 2);
      my $basename = basename("$fname",  ".faa");
      my ($base, $path, $suffix) = fileparse($fname, qr/\.\w+\.faa/);
      print {$out} "NAME  $base\nACC   KO:", join(",",  @{$ortho_rev{$base}} ), "\n";
      my @desc;
      foreach my $uniprot (@{$ortho_rev{$base}}) {
         push (@desc, $keg{$uniprot});
      }
      print {$out} "DESC  ", join(",", @desc), "\n"; 
   }
   unless ($_ =~ m/^NAME/){
   print {$out} $_;}
}
# run it either twice with and without CYCLE METAB or /
#system("cp $meta.new $meta.alt.new"); afterwards remove first line from file
#system("perl -i.bak -ne 'print if $. > 1' meta.hmm.new"); <- use this for hmmpress> File becomes to big to load nano or gedit
#to add a line to the beginning of the file "sed -i '1s/^/#CYCLE METAB\n/' file"