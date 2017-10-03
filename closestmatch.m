%3 part 2
function [human,other] = closestmatch(accession)
N=50;
orig_acc=getgenbank(char(accession));
a = topnhits(accession, N);
human=[];
other=[];
%first hit returned is the same as the original query, so ingoring this
for i=2:N
    gb_data=getgenbank(char(a(i)));
    if strfind(gb_data.Source,'Homo sapiens')
        human=a(i);
        break;
    end
end
for i=2:N
    gb_data=getgenbank(char(a(i)));
    if isempty(strfind(gb_data.Source,'Homo sapiens')) && isempty(strfind(gb_data.Source, orig_acc.Source))
        other=a(i);
        break;
    end
end
if isempty(human)
    human='No human match';
end
end
