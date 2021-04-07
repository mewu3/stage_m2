#!/usr/bin/perl -w

# example: fasta_tile_overlap.pl <input_fasta> <split_size> <overlap> <region_file>
# generates input_fasta.part1, input_fasta.part2, etc, splitting each sequence in the fasta file or fasta alignment into overlapping chunks of size split_size and overlap
# Last file is optional if you want to cover only selected regions, it will only tile across these regions. If this file is not present, it will tile across the whole genome.  Example format is below,but  w/out the #, space or tab separated region name, start position, end position. It will try to find primers within the overlap distance of the ends of the region.If region is larger than split size, it will tile multiple chunks over the region. If region is smaller than the split size, it will just use the region size.
#Nsp1	100	1100
#Nsp3	2600	8850
#Nsp14	18050	19850
#3primeGenes	21600	30550

my $infile     = $ARGV[0];
my $original_split_size= $ARGV[1];
my $overlap = $ARGV[2];
my $region_file = "";
if ($ARGV[3]) {
    $region_file = $ARGV[3];
}

my @ids=();
my %strain=&read_fasta_input();

my $out;
my $longest=0;
foreach my $id (@ids) {
    if (  length($strain{$id}) > $longest ) {
	$longest = length($strain{$id});
    }
}

#print "longest: $longest\n";
my %regions=();
if ($region_file ne "" && -s $region_file) {
    open IN,"$region_file";
    my @lines=<IN>;
    chomp @lines;
    close IN;
    foreach my $line (@lines) {
	my ($name,$start,$end)=split/\s+/,$line;
	@{$regions{$name}}=($start,$end);
    }
} else {
    my $name="all";
    @{$regions{$name}}=(0,$longest);
}

foreach my $region (keys %regions) {
    print "$region\t@{$regions{$region}}\n";
}

foreach my $region (keys %regions) {
    my $split_size=$original_split_size;
    my ($start,$end)=@{$regions{$region}};
    my $region_len = $end - $start;
    my $remainder = $region_len % ($split_size-$overlap);
    print STDERR "remainder: $remainder\n";
    my $num_parts = int( $region_len / ($split_size-$overlap) );


    if ($remainder <= $split_size/2 && $num_parts > 0) {
	my $new_split_size=int($remainder/$num_parts) + 1 + $split_size;
	print STDERR "Changing split_size from $split_size to $new_split_size so there will be no remainder.\n";
	$split_size = $new_split_size;
    } else {
	$num_parts++;  # If remainder is big, add an extra file to cover the last piece
    }

    print STDERR "region $region num_parts: $num_parts\n";

    my $N="N"x(int($overlap/2+40));
    for (my $i=0; $i < $num_parts; $i++) {

	$out = $region.".".$i."part";
	open OUT,">$out" or die "Can't open $out: $!\n";
	foreach my $id (@ids) {
	    my $begin=$start + $i*($split_size-$overlap);
	    my $finish=$begin+$split_size;
	    my $part_len;
	    if ($finish>$end) {
		$part_len=$end-$begin+1;
		$finish=$end;
	    } else {
		$part_len=$split_size;
	    }
        print "$begin $finish\n";
	    my $seq=substr($strain{$id}, $start + $i*($split_size-$overlap),$part_len);
	    $seq =~ s/[-|\*]//g;
	    $seq =~ s/[^atcgATCG]/N/g;
	    my $piece;
	    if ($overlap/2 >= $split_size) {
		$piece=$split_size/2;
	    } else {
		$piece = int($overlap/2)+1;
	    }
	    my $front=substr($seq,0,$piece);
	    my $back=substr($seq,-$piece);
	    my $glued_ends=$front.$N.$back;
	    if ($overlap > 36) {
		print OUT ">$id alignmentPositions $begin-$finish\n$glued_ends\n"; # You want primers toward the ends of the regions, so you can get overlapping amplicons from each region.
	    } else {
		print OUT ">$id alignmentPositions $begin-$finish\n$seq\n";  # the overlap is too small to get primers that don't overlap from different regions, so just assume the user simply wants to prime anywhere on the region, and regions won't necessarily overlap.
	    }
	}
	close OUT;
    } # for (my $i=0; $i < $num_parts; $i++) {
} # foreach my $region (keys %regions) {


sub read_fasta_input {

    open INDATA,"$infile" or die "Can't open $infile: $!\n";

    #get all sequences and make a hash called strain with id and sequence info
    my $sequence = "";
    my $id;
    my $num_genomes=0;
    my %strain=();
    while (my $line = <INDATA>) {
        #chomp $line;
        if ($line =~ /^>(.*)/) {
            if ($sequence ne "") {
                $ids[$num_genomes]=$id;
		$strain{$id}=uc($sequence); #make hash with $id as keys and $sequence as values
                $num_genomes++;
            }
            $id = $1;
            $id =~ s/\r//g;
            $id =~ s/\t/ /g;
            $id =~ s/\s+/ /g;
            chomp($id);
            $sequence = "";
        } else {
            chomp($line);
	    my $n10="N"x10;
           # $line =~ s/[-|\*]//g; #Get rid of "-" and "*" if it's from an alignment
            $line =~ s/[\n\t\s\r]//g; # Get rid of spaces, tabs, and new lines
            #$line =~ s/[^atcgATCG]/N/g; # replace non-ATCG chars with N
	    #$line =~ s/N{10,}/$n10/g; # replace 10 or more N's with just 10 N's
            $sequence .= $line;
        }

    }

    # Add last one
    if (defined $id) {
        $ids[$num_genomes]=$id;
        $strain{$id}=uc($sequence);
        $num_genomes++;
    }

    close INDATA or warn  "Can't close $infile: $!\n";
    #print STDERR "There are $num_genomes in fasta file.\n";
    return %strain;

} # end sub read_fasta_input
