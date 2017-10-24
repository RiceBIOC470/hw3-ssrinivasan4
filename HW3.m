% GB comments
1.	100
2a. 70 Be careful with the use of the function showalignment. Feeding your align into the function only gives a snippet of the entire possible coding sequence and therefore outputs a artificially high percent alignment.
2b. 70 same issue as 2a. 
2c. 70 same issue as 2a
3a 100
3b. 100 
3c. 100  	
Overall: 87


%HW3
%Sanjana Srinivasan

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 
matchval=2;
mismatchval=-1;
gapval=-1;
ofdiag=ones(7) - eye(7);
S=matchval*eye(7)+mismatchval*ofdiag;

seq1='GTAATCC';
seq2='GTATCCG';

figure;
[score, align, start]=swalign(seq1, seq2, 'Alphabet', 'nt', 'ScoringMatrix', S, 'GapOpen', 1, 'Showscore', true);
align;


%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 
erk1=getgenbank('NM_002746');
erk2=getgenbank('NM_002745');
coding_erk1=erk1.CDS.indices;
coding_erk2=erk2.CDS.indices;
[score_dna,align_dna]=swalign(erk1.Sequence(coding_erk1(1):coding_erk1(2)), erk2.Sequence(coding_erk2(1):coding_erk2(2)));
showalignment(align_dna);
%there is 819/1088 (75%) identity between the two sequences, and
%897/1088 (82%) positive matches between the two sequences. 


% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?
[score_aa, align_aa]=swalign(erk1.CDS.translation, erk2.CDS.translation);
showalignment(align_aa);
%there is 305/346 (88%) identity between the two amino acid sequences, and
%332/346 (96%) positive matches between the two amino acid sequences. 

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 
erk1_mouse=getgenbank('NM_011952');
merk1_coding=erk1_mouse.CDS.indices;
erk2_mouse=getgenbank('NM_011949');
merk2_coding=erk2_mouse.CDS.indices;

[erk1_score, erk1_align]=swalign(erk1.Sequence(coding_erk1(1):coding_erk1(2)), erk1_mouse.Sequence(merk1_coding(1):merk1_coding(2)));
showalignment(erk1_align);

%there is 1027/1138 (90%) identity between the two erk1 dna sequences, and
%1052/1138 (92%) positive matches. 

[erk2_score, erk2_align]=swalign(erk2.Sequence(coding_erk2(1):coding_erk2(2)), erk2_mouse.Sequence(merk2_coding(1):merk2_coding(2)));
showalignment(erk2_align);

%there is 998/1078 (93%) identity between the two erk2 dna sequences, and
%1026/1078 (95%) positive matches. 


[erk1_score_aa, erk1_align_aa]=swalign(erk1.CDS.translation, erk1_mouse.CDS.translation);
showalignment(erk1_align_aa);

%there is 367/378 (97%) identity between the two erk1 amino acid sequences, and
%370/378 (98%) positive matches.

[erk2_score_aa, erk2_align_aa]=swalign(erk2.CDS.translation, erk2_mouse.CDS.translation);
showalignment(erk2_align_aa);

%there is 355/357 (99%) identity between the two erk1 amino acid sequences, and
%357/357 (100%) positive matches.


%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 


accession='NM_002746';
N=10;

acc_num=topnhits(accession, N);

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

[human, other] = closestmatch(accession,N);

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

%human
accession='NM_033360 '; %kras variant a mrna
[human, other] = closestmatch(accession);
%the closest scoring human match was the kras mrna transcript variant b, which
%makes sense as we would expect it to be the most closely related gene with
%regards to identity
%the closest non-human match was the kras mrna transcript variant 1 in the
%gorilla, which would also make sense given that gorillas are closely
%related to humans.
accession='XR_622479'; %predicted braf gene in the common marmoset
[human, other] = closestmatch(accession);

%the closest scoring human gene is complete coding sequence of the braf
%gene. One would assume that no high scoring human gene would be found
%in the top 50 hits if the organism chosen was not a closely related
%primate. 

%the closest non-human scoring gene is the predicted braf transcript
%variant X4 in the black capped squirrel monkey. This makes sense given
%that both are primates.
