%3 part 1


function a = topnhits(accession, N)
gb_dat=getgenbank(char(accession));
[RID, ROTE] = blastncbi(gb_dat.Sequence, 'blastn');
blast_data=getblast(RID, 'WaitTime', ROTE);
for i=1:N
blast_data_hits=blast_data.Hits(i).Name;
split_hits=strsplit(blast_data_hits, '|');
a(1,i)=split_hits(4);
end
end